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
  title: [ $V=4pi\/3 * 10^3 upright(f m)^3$ ],
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
  title: [ $V=4pi\/3 * 100^3 upright(f m)^3$ ],
  xlim: (10.0, 300.0),
  ylim: (-2, 0.1),
  xaxis:(subticks:1),
  yaxis:(subticks:9, ticks:range(-2, 0, step:1)),
  legend: (position: top + left),
  ..plots_eV1
)
]

#pagebreak()





#let data_cs1 = csv("cs/rep_T=110.0_R=30.0_e=0.0_all.csv")
#let data_cs2 = csv("cs/rep_T=110.0_R=30.0_e=0.3_all.csv")
#let data_cs3 = csv("cs/rep_T=110.0_R=30.0_e=0.7_all.csv")
#let data_cs4 = csv("cs/rep_T=110.0_R=30.0_e=1.0_all.csv")

// 提取数据列
#let T_data1 = data_cs1.slice(1).map(row => float(row.at(1)))
#let cs2_data1 = data_cs1.slice(1).map(row => float(row.at(4)))

#let T_data2 = data_cs2.slice(1).map(row => float(row.at(1)))
#let cs2_data2 = data_cs2.slice(1).map(row => float(row.at(1)))

#let T_data3 = data_cs3.slice(1).map(row => float(row.at(1)))
#let cs2_data3 = data_cs3.slice(1).map(row => float(row.at(4)))

#let T_data4 = data_cs4.slice(1).map(row => float(row.at(1)))
#let cs2_data4 = data_cs4.slice(1).map(row => float(row.at(4)))


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

#pagebreak()



#let data_cs31a = csv("cs/rep_T=90.0_R=30.0_e=0.0_before.csv")
#let data_cs31b = csv("cs/rep_T=90.0_R=30.0_e=0.0_after.csv")

#let data_cs32a= csv("cs/rep_T=90.0_R=30.0_e=0.3_before.csv")
#let data_cs32b = csv("cs/rep_T=90.0_R=30.0_e=0.3_after.csv")

#let data_cs33 = csv("cs/rep_T=90.0_R=30.0_e=0.7_all.csv")
#let data_cs34 = csv("cs/rep_T=90.0_R=30.0_e=1.0_all.csv")



// 提取数据列
#let T_data31a = data_cs31a.slice(1).map(row => float(row.at(1)))
#let T_data31b = data_cs31b.slice(1).map(row => float(row.at(1)))
#let T_data32a = data_cs32a.slice(1).map(row => float(row.at(1)))
#let T_data32b = data_cs32b.slice(1).map(row => float(row.at(1)))
#let T_data33 = data_cs33.slice(1).map(row => float(row.at(1)))
#let T_data34 = data_cs34.slice(1).map(row => float(row.at(1)))



#let cs_31a = data_cs31a.slice(1).map(row => float(row.at(4)))
#let cs_31b = data_cs31b.slice(1).map(row => float(row.at(4)))

#let cs_32a = data_cs32a.slice(1).map(row => float(row.at(4)))
#let cs_32b = data_cs32b.slice(1).map(row => float(row.at(4)))

#let cs_33 = data_cs33.slice(1).map(row => float(row.at(4)))
#let cs_34 = data_cs34.slice(1).map(row => float(row.at(4)))

#let cv_31a = data_cs31a.slice(1).map(row => float(row.at(5)))
#let cv_31b = data_cs31b.slice(1).map(row => float(row.at(5)))
#let cv_32a = data_cs32a.slice(1).map(row => float(row.at(5)))
#let cv_32b = data_cs32b.slice(1).map(row => float(row.at(5)))
#let cv_33 = data_cs33.slice(1).map(row => float(row.at(5)))
#let cv_34 = data_cs34.slice(1).map(row => float(row.at(5)))







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
    legend: (position: top + right),
    lq.plot(
      T_data31a, cs_31a,
      stroke: (thickness:2pt, paint: blue, dash: lines_style.b),
      mark: none,
      label: [ $e=0.0$ ]
    ),

    lq.plot(
      T_data31b, cs_31b,
      stroke: (thickness:2pt, paint: blue, dash: lines_style.b),
      mark: none,
      //label: [ $e=0.0$ ]
    ),
    lq.plot(
      T_data32a, cs_32a,
      stroke: (thickness:2pt, paint: red, dash: lines_style.b),
      mark: none,
      label: [ $e=0.3$ ]
    ),
    lq.plot(
      T_data32b, cs_32b,
      stroke: (thickness:2pt, paint: red, dash: lines_style.b),
      mark: none,
      //label: [ $e=0.3$ ]
    ),
    lq.plot(
      T_data33, cs_33,
      stroke: (thickness:2pt, paint: green, dash: lines_style.c),
      mark: none,
      label: [ $e=0.7$ ]
    ),
    lq.plot(
      T_data34, cs_34,
      stroke: (thickness:2pt, paint: orange, dash: lines_style.d),
      mark: none,
      label: [ $e=1.0$ ]
      )
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
    ylim: (0.0, 4.0),
    xaxis:(subticks:1),
    yaxis:(subticks:1, ticks:make-yaxis(0.0, 4.0, 0.5)),
    legend: (position: top + left),
    lq.plot(
      T_data31a, cv_31a,
      stroke: (thickness:2pt, paint: blue, dash: lines_style.b),
      mark: none,
      label: [ $e=0.0$ ]
    ),

    lq.plot(
      T_data31b, cv_31b,
      stroke: (thickness:2pt, paint: blue, dash: lines_style.b),
      mark: none,
      //label: [ $e=0.0$ ]
    ),
    lq.plot(
      T_data32a, cv_32a,
      stroke: (thickness:2pt, paint: red, dash: lines_style.b),
      mark: none,
      label: [ $e=0.3$ ]
    ),
    lq.plot(
      T_data32b, cv_32b,
      stroke: (thickness:2pt, paint: red, dash: lines_style.b),
      mark: none,
      //label: [ $e=0.3$ ]
    ),
    lq.plot(
      T_data33, cv_33,
      stroke: (thickness:2pt, paint: green, dash: lines_style.c),
      mark: none,
      label: [ $e=0.7$ ]
    ),
    lq.plot(
      T_data34, cv_34,
      stroke: (thickness:2pt, paint: orange, dash: lines_style.d),
      mark: none,
      label: [ $e=1.0$ ]
      )
  )
]


#pagebreak()