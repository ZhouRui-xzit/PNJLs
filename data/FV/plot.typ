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

== 声速曲线


#let data_cs1 = csv("rep_mu=R30.0e0.0_950.0.csv")
#let data_cs2 = csv("rep_mu=R30.0e0.3_950.0.csv")
#let data_cs3 = csv("rep_mu=R30.0e0.7_950.0.csv")
#let data_cs4 = csv("rep_mu=R30.0e1.0_950.0.csv")

// 提取数据列
#let T_data1 = data_cs1.slice(1).map(row => float(row.at(0)))
#let cs2_data1 = data_cs1.slice(1).map(row => float(row.at(4)))

#let T_data2 = data_cs2.slice(1).map(row => float(row.at(0)))
#let cs2_data2 = data_cs2.slice(1).map(row => float(row.at(4)))

#let T_data3 = data_cs3.slice(1).map(row => float(row.at(0)))
#let cs2_data3 = data_cs3.slice(1).map(row => float(row.at(4)))

#let T_data4 = data_cs4.slice(1).map(row => float(row.at(0)))
#let cs2_data4 = data_cs4.slice(1).map(row => float(row.at(4)))






#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $c_(s\/rho_B)^2$ ],
    title: [  $mu_B = 950 upright(M e V)$  ],
    xlim: (50.0, 300.0),
    ylim: (0.0, 0.35),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: bottom + right),
    lq.plot(
      T_data1, cs2_data1,
      stroke: (thickness:2pt, paint: blue, dash: lines_style.b),
      mark: none,
      label: [ $e=0.0$ ]
    ),

    lq.plot(
      T_data2, cs2_data2,
      stroke: (thickness:2pt, paint: red, dash: lines_style.b),
      mark: none,
      label: [ $e=0.3$ ]
    ),
    lq.plot(
      T_data3, cs2_data3,
      stroke: (thickness:2pt, paint: green, dash: lines_style.c),
      mark: none,
      label: [ $e=0.7$ ]
    ),

    lq.plot(
      T_data4, cs2_data4,
      stroke: (thickness:2pt, paint: orange, dash: lines_style.d),
      mark: none,
      label: [ $e=1.0$ ]
    )

  )
]



#let data_cs01 = csv("rep_mu=R30.0e0.0_1.0.csv")
#let data_cs02 = csv("rep_mu=R30.0e0.3_1.0.csv")
#let data_cs03 = csv("rep_mu=R30.0e0.7_1.0.csv")
#let data_cs04 = csv("rep_mu=R30.0e1.0_1.0.csv")

// 提取数据列
#let T_data01 = data_cs01.slice(1).map(row => float(row.at(0)))
#let cs2_data01 = data_cs01.slice(1).map(row => float(row.at(4)))

#let T_data02 = data_cs02.slice(1).map(row => float(row.at(0)))
#let cs2_data02 = data_cs02.slice(1).map(row => float(row.at(4)))

#let T_data03 = data_cs03.slice(1).map(row => float(row.at(0)))
#let cs2_data03 = data_cs03.slice(1).map(row => float(row.at(4)))

#let T_data04 = data_cs04.slice(1).map(row => float(row.at(0)))
#let cs2_data04 = data_cs04.slice(1).map(row => float(row.at(4)))


#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $c_(s\/rho_B)^2$ ],
    title: [  $mu_B = 0 upright(M e V)$  ],
    xlim: (50.0, 300.0),
    ylim: (0.0, 0.35),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: bottom + right),
    lq.plot(
      T_data01, cs2_data01,
      stroke: (thickness:1pt, paint: blue, dash: lines_style.b),
      mark: none,
      label: [ $e=0.0$ ]
    ),

    lq.plot(
      T_data02, cs2_data02,
      stroke: (thickness:1pt, paint: red, dash: lines_style.c),
      mark: none,
      label: [ $e=0.3$ ]
    ),
    lq.plot(
      T_data03, cs2_data03,
      stroke: (thickness:1pt, paint: green, dash: lines_style.d),
      mark: none,
      label: [ $e=0.7$ ]
    ),

    lq.plot(
      T_data04, cs2_data04,
      stroke: (thickness:1pt, paint: orange, dash: lines_style.e),
      mark: none,
      label: [ $e=1.0$ ]
    )

  )
]