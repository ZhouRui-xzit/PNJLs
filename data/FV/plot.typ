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
#let all_R_values1a = get_all_R_values("FV/T_mu0_sph_m.csv")

#let data1a = fig_phiT("FV/T_mu0_sph_m.csv")

#let plots1a = ()
#for (i, r) in all_R_values1a.enumerate() {
  let key = str(r)
  plots1a.push(
    lq.plot(
      data1a.at(key).T, data1a.at(key).phi,
      stroke: (thickness:1.5pt, paint: line_color.at(("a", "b", "c", "d", "e", "f", "g", "h").at(calc.rem(i, 8))), 
      dash:lines_style.b),
      mark: none,
      label:[ $R=#r$ ]
    )
  )
}



// 自动获取所有R值
#let all_R_values1b = get_all_R_values("FV/T_mu0_el_m.csv")

#let data1b = fig_phiT("FV/T_mu0_el_m.csv")

#let plots1b = ()

#for (i, r) in all_R_values1b.enumerate() {
  let key = str(r)
  plots1b.push(
    lq.plot(
      data1b.at(key).T, data1b.at(key).phi,
      stroke: (thickness:1.5pt, paint: line_color.at(("a", "b", "c", "d", "e", "f", "g", "h").at(calc.rem(i, 8))), 
      dash:lines_style.c),
      mark: none,
      label:[ $cal(R)=#r$ ]
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
  title: [ 真实质量MRE ],
  xlim: (10.0, 300.0),
  ylim: (-2, 0.1),
  xaxis:(subticks:1),
  yaxis:(subticks:9, ticks:range(-2, 0, step:1)),
  legend: (position: top + right),
  ..plots1a,  // 使用展开语法将所有图线添加到图表
  ..plots1b,
)
]


// 自动获取所有R值
#let all_R_values2a = get_all_R_values("FV/T_mu0_sph_D.csv")

#let data2a = fig_phiT("FV/T_mu0_sph_D.csv")

#let plots2a = ()
#for (i, r) in all_R_values2a.enumerate() {
  let key = str(r)
  plots2a.push(
    lq.plot(
      data2a.at(key).T, data2a.at(key).phi,
      stroke: (thickness:1.5pt, paint: line_color.at(("a", "b", "c", "d", "e", "f", "g", "h").at(calc.rem(i, 8))), 
      dash:lines_style.b),
      mark: none,
      label:[ $R=#r$ ]
    )
  )
}

#let all_R_values2b = get_all_R_values("FV/T_mu0_el_D.csv")
#let data2b = fig_phiT("FV/T_mu0_el_D.csv")
#let plots2b = ()
#for (i, r) in all_R_values2b.enumerate() {
  let key = str(r)
  plots2b.push(
    lq.plot(
      data2b.at(key).T, data2b.at(key).phi,
      stroke: (thickness:1.5pt, paint: line_color.at(("a", "b", "c", "d", "e", "f", "g", "h").at(calc.rem(i, 8))), 
      dash:lines_style.c),
      mark: none,
      label:[ $cal(R)=#r$ ]
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
  title: [Dirichlet 边界条件],
  xlim: (10.0, 300.0),
  ylim: (-2, 0.1),
  xaxis:(subticks:1),
  yaxis:(subticks:9, ticks:range(-2, 0, step:1)),
  legend: (position: top + right),
  ..plots2a,  // 使用展开语法将所有图线添加到图表
  ..plots2b,
)
]

// 自动获取所有R值
#let all_R_values3a = get_all_R_values("FV/T_mu0_sph_N.csv")
#let data3a = fig_phiT("FV/T_mu0_sph_N.csv")
#let plots3a = ()
#for (i, r) in all_R_values3a.enumerate() {
  let key = str(r)
  plots3a.push(
    lq.plot(
      data3a.at(key).T, data3a.at(key).phi,
      stroke: (thickness:1.5pt, paint: line_color.at(("a", "b", "c", "d", "e", "f", "g", "h").at(calc.rem(i, 8))), 
      dash:lines_style.b),
      mark: none,
      label:[ $R=#r$ ]
    )
  )
}
#let all_R_values3b = get_all_R_values("FV/T_mu0_el_N.csv")
#let data3b = fig_phiT("FV/T_mu0_el_N.csv")
#let plots3b = ()
#for (i, r) in all_R_values3b.enumerate() {
  let key = str(r)
  plots3b.push(
    lq.plot(
      data3b.at(key).T, data3b.at(key).phi,
      stroke: (thickness:1.5pt, paint: line_color.at(("a", "b", "c", "d", "e", "f", "g", "h").at(calc.rem(i, 8))), 
      dash:lines_style.c),
      mark: none,
      label:[ $cal(R)=#r$ ]
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
  title: [Neumann 边界条件],
  xlim: (10.0, 300.0),
  ylim: (-2, 0.1),
  xaxis:(subticks:1),
  yaxis:(subticks:9, ticks:range(-2, 0, step:1)),
  legend: (position: top + right),
  ..plots3a,  // 使用展开语法将所有图线添加到图表
  ..plots3b,
)
]


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
#let all_R_values_eV2 = get_all_R_values("FV/equal_VR=7.0.csv")
#let data_eV2 = fig_phiT("FV/equal_VR=7.0.csv")
#let plots_eV2 = ()
#for (i, r) in all_R_values_eV2.enumerate() {
  let key = str(r)
  plots_eV2.push(
    lq.plot(
      data_eV2.at(key).T, data_eV2.at(key).phi,
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
  title: [ $V=4pi\/3 * 7^3 upright(f m)^3$ ],
  xlim: (10.0, 300.0),
  ylim: (-2, 0.1),
  xaxis:(subticks:1),
  yaxis:(subticks:9, ticks:range(-2, 0, step:1)),
  legend: (position: top + left),
  ..plots_eV2
)
]





// 自动获取所有R值
#let all_R_values_eV3 = get_all_R_values("FV/equal_V_R=10.0.csv")
#let data_eV3 = fig_phiT("FV/equal_V_R=10.0.csv")
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
  title: [ $V=4pi\/3 * 50^3 upright(f m)^3$ ],
  xlim: (10.0, 300.0),
  ylim: (-2, 0.1),
  xaxis:(subticks:1),
  yaxis:(subticks:9, ticks:range(-2, 0, step:1)),
  legend: (position: top + left),
  ..plots_eV3
)
]