#import "@preview/lilaq:0.5.0" as lq
#import "@preview/codly:1.3.0": * // Code

#import "@preview/physica:0.9.5": *  // 物理符号和单位包
#import "@preview/numty:0.0.5" as nt  // 数值计算包

#import  "../plot_config.typ": *

// 设置正文字体为Arial（英文）和方正书宋（中文）
#set text(font: ("Arial", "FZShuSong-Z01"), size: 12pt)

// 设置粗体字体为Arial（英文）和方正黑体（中文）
#show strong: text.with(font: ("Arial","FZHei-B01"), size: 12pt)

// 设置斜体字体为Arial（英文）和方正楷体（中文）
#show emph: text.with(font: ("Arial", "FZKai-Z03"), size: 12pt)

// 设置数学公式字体
#show math.equation: set text(font: "Latin Modern Math")

// 设置页面为自动大小
#set page(width: auto, height: auto, margin: 1em)


#let make-yaxis(ymin, ymax, step) = {
  let n = int((ymax - ymin) / step)
  range(0, n + 1, step: 1)
}




// 自动获取所有R值
#let all_R_values_eV1 = get_all_R_values("FV/1st/D_V_R=30.0.csv")
#let data_eV1 = fig_phiT("FV/1st/D_V_R=30.0.csv")
#let plots_eV1 = ()
#for (i, r) in all_R_values_eV1.enumerate() {
  let key = str(r)
  plots_eV1.push(
    lq.plot(
      data_eV1.at(key).T, data_eV1.at(key).phi,
      stroke: (thickness:1.5pt, paint: line_color.at(("a", "b", "c", "d", "e", "f", "g", "h").at(calc.rem(i, 8))), 
      dash:lines_style.b),
      mark: none,
      label:[ $e=#r$ ]
    )
  )
}


// 然后将这些绘图命令作为参数传递给 diagram
#align(center)[
  #lq.diagram(
  width: 10cm,
  height: 7cm,
  xlabel: [ $T$ (MeV) ],
  ylabel: [ $phi.alt_u$ ],
  title: [ $V=4pi\/3 times  30^3 upright(f m)^3$ ],
  xlim: (10.0, 300.0),
  ylim: (-2, 0.1),
  xaxis:(subticks:1),
  yaxis:(subticks:9, ticks:range(-2, 0, step:1)),
  legend: (position: top + left),
  ..plots_eV1
)
]

#pagebreak()


// 自动获取所有R值
#let all_R_values_eV1 = get_all_R_values("FV/1st/D_V_R=100.0.csv")
#let data_eV1 = fig_phiT("FV/1st/D_V_R=100.0.csv")
#let plots_eV1 = ()
#for (i, r) in all_R_values_eV1.enumerate() {
  let key = str(r)
  plots_eV1.push(
    lq.plot(
      data_eV1.at(key).T, data_eV1.at(key).phi,
      stroke: (thickness:1.5pt, paint: line_color.at(("a", "b", "c", "d", "e", "f", "g", "h").at(calc.rem(i, 8))), 
      dash:lines_style.b),
      mark: none,
      label:[ $e=#r$ ]
    )
  )
}


// 然后将这些绘图命令作为参数传递给 diagram
#align(center)[
  #lq.diagram(
  width: 10cm,
  height: 7cm,
  xlabel: [ $T$ (MeV) ],
  ylabel: [ $phi.alt_u$ ],
  title: [ $V=4pi\/3 times  100^3 upright(f m)^3$ ],
  xlim: (10.0, 300.0),
  ylim: (-2, 0.1),
  xaxis:(subticks:1),
  yaxis:(subticks:9, ticks:range(-2, 0, step:1)),
  legend: (position: top + left),
  ..plots_eV1
)
]

#pagebreak()

// 更好的方案：扩展 config 格式，添加 split 标记
#let load_cs_data(T, R, configs) = {
  // configs: 数组，每个元素是 (e值, 线型, 颜色, 是否分段)
  // 例如: ((0.0, "b", "a", true), (0.3, "b", "b", false), ...)
  
  let result = (:)
  
  for config in configs {
    let (e, line_style, color, has_split) = if config.len() == 4 {
      config
    } else {
      // 默认没有分段
      (config.at(0), config.at(1), config.at(2), false)
    }
    
    let base_path = "cs/rep_T=" + repr(T) + "_R=" + repr(R) + "_e=" + repr(e)
    let segments = ()
    
    if has_split {
      // 读取 _before 和 _after
      let before_data = csv(base_path + "_before.csv")
      let after_data = csv(base_path + "_after.csv")
      
      segments.push((
        mu: before_data.slice(1).map(row => float(row.at(1))),
        cs2: before_data.slice(1).map(row => float(row.at(4))),
        cv: before_data.slice(1).map(row => float(row.at(5)))
      ))
      
      segments.push((
        mu: after_data.slice(1).map(row => float(row.at(1))),
        cs2: after_data.slice(1).map(row => float(row.at(4))),
        cv: after_data.slice(1).map(row => float(row.at(5)))
      ))
    } else {
      // 读取 _all
      let data = csv(base_path + "_all.csv")
      
      segments.push((
        mu: data.slice(1).map(row => float(row.at(1))),
        cs2: data.slice(1).map(row => float(row.at(4))),
        cv: data.slice(1).map(row => float(row.at(5)))
      ))
    }
    
    result.insert(str(e), (
      segments: segments,
      style: lines_style.at(line_style),
      color: line_color.at(color),
      e: e
    ))
  }
  
  return result
}

// 方法1：手动指定颜色和线型
#let cs_configs = (
  (0.0, "b", "a"),   // e=0.0, 虚线, 蓝色
  (0.3, "b", "b"),   // e=0.3, 虚线, 红色
  (0.7, "c", "c"),   // e=0.7, 点线, 绿色
  (1.0, "d", "d"),   // e=1.0, 点划线, 橙色
)

#let cs_configs2 = (
  (0.0, "b", "a", true),   // e=0.0, 虚线, 蓝色
  (0.3, "b", "b", true),   // e=0.3, 虚线, 红色
  (0.7, "c", "c"),   // e=0.7, 点线, 绿色
  (1.0, "d", "d"),   // e=1.0, 点划线, 橙色
)


#let cs_configs4 = (
  (0.0, "b", "a", true),   // e=0.0, 虚线, 蓝色
  (0.3, "b", "b", true),   // e=0.3, 虚线, 红色
  (0.7, "c", "c", true),   // e=0.7, 点线, 绿色
  (1.0, "d", "d", true),   // e=1.0, 点划线, 橙色
)

#let cs_data1 = load_cs_data(50.0, 30.0, cs_configs4)
#let cs_data2 = load_cs_data(90.0, 30.0, cs_configs2)

#let cs_data3 = load_cs_data(110.0, 30.0, cs_configs)
#let cs_data4 = load_cs_data(150.0, 30.0, cs_configs)


#let cs_data21 = load_cs_data(50.0, 100.0, cs_configs4)
#let cs_data22 = load_cs_data(123.0, 100.0, cs_configs2)

#let cs_data23 = load_cs_data(150.0, 100.0, cs_configs)





#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $mu_B$ (MeV) ],
    ylabel: [ $c_(s\/rho_B)^2$ ],
    title: [  $T = 50 upright(M e V)$  ],
    xlim: (100.0, 1200.0),
    ylim: (0.0, 0.35),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data1.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.mu, seg.cs2,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style),
          mark: none,
          label: if i == 0 { [ $e=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )
]
#pagebreak()


#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $mu_B$ (MeV) ],
    ylabel: [ $c_(s\/rho_B)^2$ ],
    title: [  $T = 90 upright(M e V)$  ],
    xlim: (100.0, 1200.0),
    ylim: (0.0, 0.35),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data2.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.mu, seg.cs2,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style),
          mark: none,
          label: if i == 0 { [ $e=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )
]
#pagebreak()

#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $mu_B$ (MeV) ],
    ylabel: [ $c_(s\/rho_B)^2$ ],
    title: [  $T = 110 upright(M e V)$  ],
    xlim: (100.0, 1200.0),
    ylim: (0.0, 0.35),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + right),
    
    // 遍历字典，每个 e 值可能有多段
    ..cs_data3.values().map(d => {
      // 为每个 segment 生成一个 plot，但只有第一个有 label
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.mu, seg.cs2,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style),
          mark: none,
          label: if i == 0 { [ $e=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )
]


#pagebreak()
#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $mu_B$ (MeV) ],
    ylabel: [ $c_(s\/rho_B)^2$ ],
    title: [  $T = 150 upright(M e V)$  ],
    xlim: (100.0, 1200.0),
    ylim: (0.0, 0.35),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + right),
    
    // 遍历字典，每个 e 值可能有多段
    ..cs_data4.values().map(d => {
      // 为每个 segment 生成一个 plot，但只有第一个有 label
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.mu, seg.cs2,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style),
          mark: none,
          label: if i == 0 { [ $e=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )
]
#pagebreak()

#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $mu_B$ (MeV) ],
    ylabel: [ $C_V$ ],
    title: [  $T = 50 upright(M e V)$  ],
    xlim: (100.0, 1200.0),
    ylim: (0.0, 1.0),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data1.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.mu, seg.cv,
          stroke: (thickness: 1.1pt, paint: d.color, dash: d.style),
          mark: none,
          label: if i == 0 { [ $e=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )
]
#pagebreak()


#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $mu_B$ (MeV) ],
    ylabel: [ $C_V$ ],
    title: [  $T = 90 upright(M e V)$  ],
    xlim: (100.0, 1200.0),
    ylim: (0.0, 3.5),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data2.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.mu, seg.cv,
          stroke: (thickness: 1.1pt, paint: d.color, dash: d.style),
          mark: none,
          label: if i == 0 { [ $e=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )
]
#pagebreak()


#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $mu_B$ (MeV) ],
    ylabel: [ $C_V$ ],
    title: [  $T = 110 upright(M e V)$  ],
    xlim: (100.0, 1200.0),
    ylim: (0.0, 6.0),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data3.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.mu, seg.cv,
          stroke: (thickness: 1.1pt, paint: d.color, dash: d.style),
          mark: none,
          label: if i == 0 { [ $e=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )
]
#pagebreak()

#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $mu_B$ (MeV) ],
    ylabel: [ $C_V$ ],
    title: [  $T = 150 upright(M e V)$  ],
    xlim: (100.0, 1200.0),
    ylim: (0.0, 20.0),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data4.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.mu, seg.cv,
          stroke: (thickness: 1.1pt, paint: d.color, dash: d.style),
          mark: none,
          label: if i == 0 { [ $e=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )
]
#pagebreak()




#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $mu_B$ (MeV) ],
    ylabel: [ $c_(s\/rho_B)^2$ ],
    title: [  $T = 50 upright(M e V)$  ],
    xlim: (100.0, 1200.0),
    ylim: (0.0, 0.35),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data21.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.mu, seg.cs2,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style),
          mark: none,
          label: if i == 0 { [ $e=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )
]
#pagebreak()

#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $mu_B$ (MeV) ],
    ylabel: [ $c_(s\/rho_B)^2$ ],
    title: [  $T = 123 upright(M e V)$  ],
    xlim: (100.0, 1200.0),
    ylim: (0.0, 0.35),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data22.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.mu, seg.cs2,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style),
          mark: none,
          label: if i == 0 { [ $e=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )
]
#pagebreak()


#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $mu_B$ (MeV) ],
    ylabel: [ $c_(s\/rho_B)^2$ ],
    title: [  $T = 150 upright(M e V)$  ],
    xlim: (100.0, 1200.0),
    ylim: (0.0, 0.35),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data23.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.mu, seg.cs2,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style),
          mark: none,
          label: if i == 0 { [ $e=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )
]
#pagebreak()



#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $mu_B$ (MeV) ],
    ylabel: [ $C_V$ ],
    title: [  $T = 50 upright(M e V)$  ],
    xlim: (100.0, 1200.0),
    ylim: (0.0, 1.0),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data21.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.mu, seg.cv,
          stroke: (thickness: 1.1pt, paint: d.color, dash: d.style),
          mark: none,
          label: if i == 0 { [ $e=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )
]
#pagebreak()


#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $mu_B$ (MeV) ],
    ylabel: [ $C_V$ ],
    title: [  $T = 123 upright(M e V)$  ],
    xlim: (100.0, 1200.0),
    ylim: (0.0, 10.0),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data22.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.mu, seg.cv,
          stroke: (thickness: 1.1pt, paint: d.color, dash: d.style),
          mark: none,
          label: if i == 0 { [ $e=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )
]
#pagebreak()



#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $mu_B$ (MeV) ],
    ylabel: [ $C_V$ ],
    title: [  $T = 150 upright(M e V)$  ],
    xlim: (100.0, 1200.0),
    ylim: (0.0, 21.0),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data23.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.mu, seg.cv,
          stroke: (thickness: 1.1pt, paint: d.color, dash: d.style),
          mark: none,
          label: if i == 0 { [ $e=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )
]
