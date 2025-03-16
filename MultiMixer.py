import re
import math
import argparse
from itertools import zip_longest

def parse_filter_conditions(filter_str):
    """解析过滤条件字符串为结构化数据"""
    conditions = []
    if not filter_str:
        return conditions
    
    for cond in filter_str.split(','):
        parts = cond.split(':')
        if len(parts) != 2:
            raise ValueError(f"无效过滤条件格式: {cond}")
        
        pdg = int(parts[0])
        expr = parts[1]
        
        # 解析比较运算符
        operators = ['>=', '<=', '>', '<']
        op = None
        for o in operators:
            if o in expr:
                op = o
                param, value = expr.split(o, 1)
                break
        
        if not op:
            raise ValueError(f"无效运算符: {expr}")
        
        conditions.append({
            'pdg': pdg,
            'param': param.strip(),
            'op': op,
            'value': float(value)
        })
    return conditions

def calculate_kinematics(px, py, pz, energy):
    """计算粒子运动学参数"""
    pT = math.sqrt(px**2 + py**2)
    try:
        eta = 0.5 * math.log((energy + pz) / (energy - pz)) if energy != pz else 0.0
    except:
        eta = 0.0
    phi = math.atan2(py, px)
    return {
        'pT': pT,
        'eta': eta,
        'phi': phi,
        'energy': energy,
        'pz': pz
    }

def check_event(event_str, conditions):
    """检查事件是否满足所有过滤条件"""
    if not conditions:
        return True
    
    event_lines = event_str.strip().split('\n')
    particles = event_lines[2:-1]  # 跳过事件头和尾
    
    for cond in conditions:
        pdg = cond['pdg']
        param = cond['param']
        op = cond['op']
        threshold = cond['value']
        
        # 查找匹配的粒子
        matched = False
        for p_line in particles:
            parts = p_line.strip().split()
            if len(parts) < 11:
                continue
            
            if int(parts[0]) != pdg:
                continue
                
            # 解析四动量
            px = float(parts[6])
            py = float(parts[7])
            pz = float(parts[8])
            energy = float(parts[9])
            
            # 计算运动学参数
            kinematics = calculate_kinematics(px, py, pz, energy)
            value = kinematics.get(param, None)
            if value is None:
                continue
                
            # 进行比较
            if op == '>' and not (value > threshold):
                return False
            elif op == '<' and not (value < threshold):
                return False
            elif op == '>=' and not (value >= threshold):
                return False
            elif op == '<=' and not (value <= threshold):
                return False
        
        # 如果该条件类型需要至少一个粒子存在
        if not any(int(p.split()[0]) == pdg for p in particles):
            return False  # 如果要求必须存在该粒子
            
    return True

def read_events(filename, max_events=None, filter_conditions=None):
    """从文件读取事件并过滤"""
    pre_lines = []
    events = []
    post_lines = []
    current_event = []
    in_event = False
    count = 0

    with open(filename, 'r') as f:
        for line in f:
            if line.strip().startswith('<event>'):
                in_event = True
                current_event = [line]
            elif line.strip().startswith('</event>'):
                current_event.append(line)
                event_str = ''.join(current_event)
                
                # 应用过滤条件
                if check_event(event_str, filter_conditions):
                    events.append(event_str)
                    count += 1
                    if max_events and count >= max_events:
                        break
                
                current_event = []
                in_event = False
            elif in_event:
                current_event.append(line)
            else:
                if not events and not in_event:
                    pre_lines.append(line)
                else:
                    post_lines.append(line)
    return {
        'pre': pre_lines,
        'events': events,
        'post': post_lines,
        'filename': filename
    }

def adjust_particles(particles, base_index, color_offset):
    """调整粒子索引和颜色标签"""
    adjusted = []
    for p in particles:
        parts = p.strip().split()
        
        # 调整母粒子索引
        mothers = []
        for idx in [2,3]:
            original = int(parts[idx])
            if original > 0:
                mothers.append(str(original + base_index))
            else:
                mothers.append(str(original))
        
        # 调整颜色标签
        colors = []
        for idx in [4,5]:
            original = parts[idx]
            if original != '0':
                colors.append(str(int(original) + color_offset))
            else:
                colors.append(original)
        
        # 重组粒子行
        new_p = parts[:2] + mothers + colors + parts[6:]
        adjusted.append('  '.join(new_p))
    return adjusted

def merge_event_group(event_group):
    """合并一个事件组（多个事件）为一个复合事件"""
    total_particles = 0
    total_weight = 1.0
    color_offset = 0
    all_particles = []
    event_headers = []

    # 收集所有事件头信息
    for event_str in event_group:
        event_lines = event_str.strip().split('\n')
        header = event_lines[1].strip().split()
        event_headers.append({
            'npart': int(header[0]),
            'weight': float(header[2]),
            'scale': header[3],
            'alpha_qed': header[4],
            'alpha_qcd': header[5]
        })

    # 处理每个事件
    for idx, event_str in enumerate(event_group):
        event_lines = event_str.strip().split('\n')
        particles = event_lines[2:-1]
        
        # 调整当前事件的粒子
        adjusted = adjust_particles(
            particles,
            base_index=total_particles,
            color_offset=color_offset
        )
        all_particles.extend(adjusted)
        
        # 更新累积参数
        total_particles += event_headers[idx]['npart']
        total_weight *= event_headers[idx]['weight']
        
        # 更新颜色偏移
        max_color = max(int(tag) for p in particles 
                      for tag in p.strip().split()[4:6] 
                      if tag != '0') or 0
        color_offset += max_color + 1

    # 构建新事件头
    new_header = [
        str(total_particles),        # 总粒子数
        '999',                        # 过程代码
        f"{total_weight:.6E}",        # 合并权重
        event_headers[0]['scale'],    # 标度
        event_headers[0]['alpha_qed'],# alpha_QED
        event_headers[0]['alpha_qcd'] # alpha_QCD
    ]

    return '\n'.join([
        '<event>',
        '  ' + '  '.join(new_header),
        *['  '.join(p.split()) for p in all_particles],
        '</event>'
    ])

def batch_generator(event_pools, batch_sizes):
    """生成事件批次的生成器"""
    iterators = []
    for pool, size in zip(event_pools, batch_sizes):
        # 将事件池按指定大小分块
        chunks = [pool[i:i+size] for i in range(0, len(pool), size)]
        iterators.append(iter(chunks))
    
    while True:
        batch = []
        for it in iterators:
            try:
                batch.extend(next(it))
            except StopIteration:
                return
        yield batch

def main(files_config, output_file, max_total_events, filter_conditions):
    """主处理函数"""
    event_pools = []
    batch_sizes = []
    
    # 计算最大可生成组数
    max_groups = min(
        len(file_info['events']) // fc['events_per_group']
        for fc, file_info in zip(files_config, [read_events(fc['filename'], fc['max_events'], filter_conditions) for fc in files_config])
    )
    
    # 加载事件池
    for fc in files_config:
        file_info = read_events(fc['filename'], fc['max_events'], filter_conditions)
        required = fc['events_per_group'] * max_groups
        event_pools.append(file_info['events'][:required])
        batch_sizes.append(fc['events_per_group'])
    
    # 生成合并事件
    generator = batch_generator(event_pools, batch_sizes)
    merged_events = []
    
    for event_group in generator:
        if len(merged_events) >= max_total_events:
            break
        merged_events.append(merge_event_group(event_group))
    
    # 写入文件
    first_file_info = read_events(files_config[0]['filename'], 0)
    with open(output_file, 'w') as f:
        f.writelines(first_file_info['pre'])
        f.write('\n'.join(merged_events))
        f.write('\n')
        f.writelines(first_file_info['post'])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='高级LHE事件混合工具')
    parser.add_argument('-o', '--output', required=True, help='输出文件路径')
    parser.add_argument('-f', '--files', nargs='+', required=True,
                       help="输入文件配置，格式为 filename:每组数量:最大事件数")
    parser.add_argument('--max-events', type=int, default=None,
                       help='最大生成事件总数')
    parser.add_argument('--filter', type=str, default=None,
                       help='粒子过滤条件，格式如 443:pT>10,13:eta<2.5')
    
    args = parser.parse_args()
    
    # 解析文件配置
    files_config = []
    for item in args.files:
        parts = item.split(':')
        if len(parts) != 3:
            raise ValueError("参数格式错误，应使用 filename:每组数量:最大事件数")
        files_config.append({
            'filename': parts[0],
            'events_per_group': int(parts[1]),
            'max_events': int(parts[2])
        })
    
    # 解析过滤条件
    filter_conditions = parse_filter_conditions(args.filter) if args.filter else []
    
    main(
        files_config=files_config,
        output_file=args.output,
        max_total_events=args.max_events if args.max_events else float('inf'),
        filter_conditions=filter_conditions
    )