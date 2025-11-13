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






= phi-T 曲线






// 自动获取所有R值
#let all_R_values_eV1 = get_all_R_values("FV/equal_VR=10.0.csv")
#let data_eV1 = fig_phiT("FV/equal_VR=10.0.csv")
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
  title: [ $V=4pi\/3 * 10^3 upright(f m)^3$ ],
  xlim: (10.0, 300.0),
  ylim: (-2, 0.1),
  xaxis:(subticks:1),
  yaxis:(subticks:9, ticks:range(-2, 0, step:1)),
  legend: (position: top + left),
  ..plots_eV1
)
]






// 自动获取所有R值
#let all_R_values_eV3 = get_all_R_values("FV/D_V_R=30.0.csv")
#let data_eV3 = fig_phiT("FV/D_V_R=30.0.csv")
#let plots_eV3 = ()
#for (i, r) in all_R_values_eV3.enumerate() {
  let key = str(r)
  plots_eV3.push(
    lq.plot(
      data_eV3.at(key).T, data_eV3.at(key).phi,
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
  title: [ $V=4pi\/3 * 30^3 upright(f m)^3$ ],
  xlim: (10.0, 300.0),
  ylim: (-2, 0.1),
  xaxis:(subticks:1),
  yaxis:(subticks:9, ticks:range(-2, 0, step:1)),
  legend: (position: top + left),
  ..plots_eV3
)
]


// 自动获取所有R值
#let all_R_values_eV4 = get_all_R_values("FV/N_V_R=10.0.csv")
#let data_eV4 = fig_phiT("FV/N_V_R=10.0.csv")
#let plots_eV4 = ()
#for (i, r) in all_R_values_eV4.enumerate() {
  let key = str(r)
  plots_eV4.push(
    lq.plot(
      data_eV4.at(key).T, data_eV4.at(key).phi,
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
  title: [ $V=4pi\/3 * 10^3 upright(f m)^3$ ],
  xlim: (10.0, 300.0),
  ylim: (-2, 0.1),
  xaxis:(subticks:1),
  yaxis:(subticks:9, ticks:range(-2, 0, step:1)),
  legend: (position: top + left),
  ..plots_eV4
)
]