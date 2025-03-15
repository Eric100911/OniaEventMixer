import re
import argparse
from itertools import zip_longest

def read_events(filename, max_events=None):
    """从指定文件读取最多max_events个事件"""
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
                events.append(''.join(current_event))
                current_event = []
                in_event = False
                count += 1
                if max_events and count >= max_events:
                    break
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

def main(files_config, output_file):
    """主处理函数"""
    # 从各文件读取事件池
    event_pools = []
    batch_sizes = []
    # 确定最终生成的合并后事件数量，考虑到每个文件的最大事件分组数，取最小值
    max_groups = min(fc['max_events'] // fc['events_per_group'] for fc in files_config)
    for fc in files_config:
        file_info = read_events(fc['filename'], fc['max_events'])
        event_pool = file_info['events'][:fc['events_per_group'] * max_groups]
        event_pools.append(event_pool)
        batch_sizes.append(fc['events_per_group'])
    
    # 创建批次生成器
    generator = batch_generator(event_pools, batch_sizes)
    
    # 合并所有批次
    merged_events = []
    for event_group in generator:
        merged_events.append(merge_event_group(event_group))
    
    # 写入文件（使用第一个文件的头信息）
    first_file_info = read_events(files_config[0]['filename'], 0)
    with open(output_file, 'w') as f:
        f.writelines(first_file_info['pre'])
        f.write('\n'.join(merged_events))
        f.write('\n')
        f.writelines(first_file_info['post'])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='LHE事件混合工具')
    parser.add_argument('-o', '--output', required=True, help='输出文件路径')
    parser.add_argument('-f', '--files', nargs='+', required=True,
                        help="输入文件配置，格式为 filename:每组数量:最大总数 如 file1.lhe:2:10 file2.lhe:1:5")
    
    args = parser.parse_args()
    
    # 解析增强型文件配置
    files_config = []
    for item in args.files:
        parts = item.split(':')
        if len(parts) != 3:
            raise ValueError("文件参数格式错误，应使用 filename:每组数量:最大总数")
        
        files_config.append({
            'filename': parts[0],
            'events_per_group': int(parts[1]),
            'max_events': int(parts[2])
        })
    
    main(files_config, args.output)