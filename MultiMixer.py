#!/usr/bin/env python3
import argparse
import json
import random
import itertools as it
from typing import Dict, List, Generator
import pylhe
from pylhe import LHEEvent, LHEEventInfo, LHEParticle, LHEFile, LHEInit
import vector

class ColorManager:
    def __init__(self):
        self._base_per_sector = 100  # 每个sector的基础间隔
        self._current_sector = -1       # 当前sector索引
        self._sector_maps = []          # 每个sector的颜色映射表

    def next_sector(self) -> None:
        """切换到下一个颜色区间"""
        self._current_sector += 1
        self._sector_maps.append({})     # 为新sector初始化空映射表

    def get_color(self, original: int) -> int:
        """获取映射后的颜色编号"""
        if original == 0:
            return 0
        # 获取当前sector的映射表
        current_map = self._sector_maps[self._current_sector]
        if original not in current_map:
            # 计算新颜色：基础值 + 原始计数器的偏移（确保每个sector独立递增）
            base = self._current_sector * self._base_per_sector
            new_color = base + 500 + len(current_map) + 1
            current_map[original] = new_color
        return current_map[original]

class ParticleSource:
    def __init__(self, config: Dict, shuffle: bool, seed: int):
        self.files = config["files"]
        self.count = config["count"]
        self.shuffle = shuffle
        self.rng = random.Random(seed)
        if shuffle:
            self.rng.shuffle(self.files)
        self.file_iter = iter(self.files)

class EventMixer:
    def __init__(self, config: Dict):
        """
        参数优先级：
        1. 显式传入的参数 (最高)
        2. JSON配置参数
        3. 默认值
        """
        # Fundamental parameters
        self.config = config['sources']
        self.max_events = config.get('max_events', 50000)
        self.filters = config.get('filters', [])
        self.shuffle = config.get('shuffle', False)
        self.seed = config.get('seed', 42)
        self.squash = config.get('squash', False)
        self.buffer = config.get('buffer', 1000)
        
        # 初始化随机数
        self.rng = random.Random(self.seed)
        self.sources = {
            name: ParticleSource(cfg, config['shuffle'], config['seed'])
            for name, cfg in config["sources"].items()
        }
        # Read init info from the first file of the first source
        # Should be compatible with one-source case
        self.lhe_file_init = pylhe.read_lhe_init((self.sources[list(self.sources.keys())[0]].files[0]))

    def to_lhe_frame(self) -> LHEFile:
        """Generate an LHEFile with anything but events"""
        lhe_file = LHEFile()
        
        # 设置初始化信息
        tmp_raw_init = self.lhe_file_init.copy()
        tmp_raw_init["nevents"] = 0
        tmp_raw_init["weightgroup"] = {} # Just to avoid pylhe error
        lhe_file.init = LHEInit(**tmp_raw_init)
        
        # 添加文件头元数据
        lhe_file.header_blocks = [
            '<!-- Merged by EventMixer -->',
            '<!-- Source config: {} -->'.format(json.dumps(self.config))
        ]
        
        return lhe_file


    def to_lhe_file(self) -> LHEFile:
        """生成完整的LHE文件对象"""
        lhe_file = LHEFile()
        
        # 设置初始化信息
        tmp_raw_init = self.lhe_file_init.copy()
        tmp_raw_init["nevents"] = self._count_total_events()
        tmp_raw_init["weightgroup"] = {} # Just to avoid pylhe error
        lhe_file.init = LHEInit(**tmp_raw_init)
        
        # 添加事件数据
        lhe_file.events = list(it.islice(self.generate(), self.max_events))
        
        # 添加文件头元数据
        lhe_file.header_blocks = [
            '<!-- Merged by EventMixer -->',
            '<!-- Source config: {} -->'.format(json.dumps(self.config))
        ]
        
        return lhe_file
    
    def _count_total_events(self) -> int:
        """预计算总事件数用于元数据"""
        try:
            return max(sum(1 for _ in self.generate()), self.max_events or 0)
        except Exception as e:
            print(f"WARNING: Event counting failed: {str(e)}")
            return 0


    def _pass_filters(self, event: LHEEvent) -> bool:
        for p in event.particles:
            for condition in self.filters:
                pid, expr = condition.split(':', 1)
                if p.id == int(pid):
                    p4 = vector.obj(px=p.px, py=p.py, pz=p.pz, e=p.e)
                    pT = p4.pt
                    eta = p4.eta
                    if not eval(expr, {'pT': pT, 'eta': eta, 'abs': abs}):
                        return False
        return True

    def _event_generator(self, files: List[str]) -> Generator[LHEEvent, None, None]:
        # Create a buffer zone.
        buffer = []
        for file in files:
            # Using itertools to extract events from the file and filter.
            buffer.extend(it.filterfalse(lambda x: not self._pass_filters(x), pylhe.read_lhe_file(file).events))
            # Once the buffer is full, do shuffling if needed, and then yield.
            if len(buffer) >= self.buffer:
                if self.shuffle:
                    self.rng.shuffle(buffer)
                yield from buffer
                buffer.clear()
            
    def _calculate_beam_total_p4(self, events: List[LHEEvent]) -> vector.MomentumObject4D:
        total_px     = sum(p.px for event in events for p in event.particles if p.status == -1)
        total_py     = sum(p.py for event in events for p in event.particles if p.status == -1)
        total_pz     = sum(p.pz for event in events for p in event.particles if p.status == -1)
        total_energy = sum(p.e  for event in events for p in event.particles if p.status == -1)
        total_p4 = vector.obj(px=total_px, py=total_py, pz=total_pz, e=total_energy)
        return total_p4

    def _squash_initial_gluons(self, events: List[LHEEvent]) -> LHEEvent:
        # 创建新的初态胶子对

        # Reassign the gluon momenta. Begin by calculating the total momentum of the beam.
        total_p4 = self._calculate_beam_total_p4(events)
        total_p3_unit = total_p4.to_3D().unit()

        # Project to the direction of p4 sum, and the redistribute the momentum.
        # Assign (p_mag - E) / 2 * p_unit and (p_mag + E) / 2 * p_unit to the two gluons as the momenta.
        gluon1_p4 = total_p3_unit.scale((total_p4.mag - total_p4.e) / 2).to_Vector4D(tau=0)
        gluon2_p4 = total_p3_unit.scale((total_p4.mag + total_p4.e) / 2).to_Vector4D(tau=0)
        

        gluon1 = LHEParticle(
            id=21, status=-1,
            mother1=0, mother2=0,
            color1=101, color2=102,
            px=gluon1_p4.px,
            py=gluon1_p4.py,
            pz=gluon1_p4.pz,
            e=gluon1_p4.e,
            m=0.0, lifetime=0.0, spin=0.0
        )
        gluon2 = LHEParticle(
            id=21, status=-1,
            mother1=0, mother2=0,
            color1=102, color2=101,
            px=gluon2_p4.px,
            py=gluon2_p4.py,
            pz=gluon2_p4.pz,
            e=gluon2_p4.e,
            m=0.0, lifetime=0.0, spin=0.0
        )

        # Dealing with colored final state particles.
        # Here's the trick: pair as much as I can, the last one connects to gluon1 and gluon2
        merged = LHEEvent(events[0].eventinfo, [])
        next_color1 = 103
        next_color2 = 104

        for event in events:
            for p in event.particles:
                if p.status == 1:
                    new_p = p
                    new_p.mother1 = 1
                    new_p.mother2 = 2
                    if p.color1 != 0 and p.color2 != 0:
                        new_p.color1 = next_color1
                        new_p.color2 = next_color2
                        if next_color1 < next_color2:
                            next_color1, next_color2 = next_color2, next_color1
                        else:
                            next_color1 += 1
                            next_color2 += 3
                    
                    merged.particles.append(new_p)

        # Unpaired coloured particles at the end of merged.particles
        if next_color1 > next_color2:
            gluon1.color2 = next_color1
            gluon2.color1 = next_color2

        merged.particles.insert(0, gluon1)
        merged.particles.insert(1, gluon2)

        # Update the number of particles
        merged.eventinfo.nparticles = len(merged.particles)

        return merged

    def _merge_events(self, events: List[LHEEvent]) -> LHEEvent:
        """处理母子关系合并的核心方法"""
        all_particles = []
        mother_index_map = {}
        color_mgr = ColorManager()

        # 第一阶段：收集所有粒子并建立索引映射
        global_index = 0
        for event_idx, event in enumerate(events):
            index_offset = global_index
            local_mother_map = {}

            # 更新颜色区间
            color_mgr.next_sector()

            # 处理每个粒子的新索引
            for part_idx, particle in enumerate(event.particles):
                new_index = global_index
                local_mother_map[part_idx + 1] = new_index + 1  # LHE索引从1开始
                global_index += 1

                # 创建新粒子（暂不处理母子关系）
                new_particle = particle
                new_particle.color1 = color_mgr.get_color(particle.color1)
                new_particle.color2 = color_mgr.get_color(particle.color2)
                all_particles.append(new_particle)
                mother_index_map[(event_idx, part_idx + 1)] = new_index + 1

            # 更新母子关系映射
            for part in all_particles[index_offset:]:
                if part.mother1 > 0:
                    part.mother1 = local_mother_map.get(part.mother1, 0)
                if part.mother2 > 0:
                    part.mother2 = local_mother_map.get(part.mother2, 0)

        # 第二阶段：跨事件母子关系处理
        for particle in all_particles:
            if particle.mother1 > 0:
                original_event = next((e_idx for e_idx, e in enumerate(events) 
                                     if particle in e.particles), 0)
                original_index = particle.mother1
                particle.mother1 = mother_index_map.get((original_event, original_index), 0)
            
            if particle.mother2 > 0:
                original_event = next((e_idx for e_idx, e in enumerate(events) 
                                     if particle in e.particles), 0)
                original_index = particle.mother2
                particle.mother2 = mother_index_map.get((original_event, original_index), 0)

        # 构建事件元数据
        base_info = events[0].eventinfo if events else LHEEventInfo(
            nparticles=0, pid=-1, weight=1.0, scale=0.0, aqed=0.007297, aqcd=0.118)
        
        event_info = LHEEventInfo(
            nparticles=len(all_particles),
            pid=base_info.pid,
            weight=sum(e.eventinfo.weight for e in events),
            scale=base_info.scale,
            aqed=base_info.aqed,
            aqcd=base_info.aqcd
        )

        return LHEEvent(eventinfo=event_info, particles=all_particles)

    def generate(self) -> Generator[LHEEvent, None, None]:
        generators = {
            name: self._event_generator(source.files)
            for name, source in self.sources.items()
        }

        # 记录活跃的生成器
        active_generators = set(generators.keys())

        while active_generators:
            batch = []
            for name in list(active_generators):  # 避免迭代时修改集合
                source = self.sources[name]
                try:
                    # 显式调用next()触发StopIteration
                    events = [next(generators[name]) for _ in range(source.count)]
                    batch.extend(events)
                except StopIteration:
                    # 生成器耗尽后移出活跃列表
                    active_generators.remove(name)

            if batch:
                if self.squash:
                    yield self._squash_initial_gluons(batch)
                else:
                    yield self._merge_events(batch)
            else:
                break  # 所有生成器已耗尽

def merge_args(config: Dict, args: argparse.Namespace) -> Dict:
    """合并命令行参数到配置字典"""
    merged = config.copy()

    if args.output is not None:
        merged['output'] = args.output
    if args.max_events is not None:
        merged['max_events'] = args.max_events
    if args.filters is not None:
        merged['filters'] = args.filters
    if args.shuffle is not None:
        merged['shuffle'] = args.shuffle
    if args.seed is not None:
        merged['seed'] = args.seed
    if args.squash is not None:
        merged['squash'] = args.squash
    if args.buffer is not None:
        merged['buffer'] = args.buffer
    
    return merged

def main():
    parser = argparse.ArgumentParser(description='LHE事件混合器')
    parser.add_argument('--config', required=True, type=str, 
                        help='JSON配置文件路径')
    parser.add_argument('--output', type=str, default=None,
                        help='覆盖配置文件中的输出路径')
    parser.add_argument('--max-events', type=int, default=None,
                        help='最大输出事件数（覆盖配置）')
    parser.add_argument('--filters', action='append', default=None,
                        help='粒子过滤条件，格式: PID:EXPR (可多次使用)')
    parser.add_argument('--shuffle', action=argparse.BooleanOptionalAction, 
                        default=None, help='是否随机混合')
    parser.add_argument('--seed', type=int, default=None,
                        help='随机数种子')
    parser.add_argument('--squash', action=argparse.BooleanOptionalAction,
                        default=None, help='压缩重复事件')
    parser.add_argument('--buffer', type=int, default=1000,
                        help='事件缓冲区大小')
    
    args = parser.parse_args()
    
    # 加载配置文件
    try:
        with open(args.config, 'r') as f:
            config = json.load(f)
    except FileNotFoundError:
        raise ValueError(f"配置文件 {args.config} 不存在")
    except json.JSONDecodeError:
        raise ValueError(f"配置文件 {args.config} 格式错误")

    # 参数合并
    merged_config = merge_args(config, args)

    output_path = merged_config.get('output', 'output.lhe')
    
    # 初始化混合器
    mixer = EventMixer(
        config=merged_config
    )

    output_file = mixer.to_lhe_file()
    pylhe.write_lhe_file_path(
        lhefile=output_file,
        filepath=output_path
    )


if __name__ == '__main__':
    main()