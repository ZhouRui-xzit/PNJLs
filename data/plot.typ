#import "@preview/lilaq:0.5.0" as lq
#import "@preview/codly:1.3.0": * // Code

#import "@preview/physica:0.9.5": *  // 物理符号和单位包
#import "@preview/numty:0.0.5" as nt  // 数值计算包

#import  "plot_config.typ": *

#set text(font: ("Libertinus Serif", "FZShuSong-Z01"), size: 12pt)
#show strong: text.with(font: ("Libertinus Serif", "FZHei-B01"), size: 12pt)
#show emph: text.with(font: ("Libertinus Serif", "FZKai-Z03"), size: 12pt)
#show math.equation: set text(font: ("Libertinus Math",), size: 12pt)

// 设置数学公式字体
#show math.equation: set text(font: "Latin Modern Math")






= M-R relation plots 
// 生成多个R-M图像的绘制命令
// 函数：批量读取CSV文件并生成绘图命令
#let generate_plots(filenames, labels) = {
  let plots = ()
  
  for (i, filename) in filenames.enumerate() {
    let data = csv(filename)
    let data_values = data.slice(1) // 跳过第一行(列名)
    
    plots.push(
      lq.plot(
        data_values.map(row => float(row.at(0))), // 第一列作为 R
        data_values.map(row => float(row.at(1))), // 第二列作为 M
        stroke: (
          thickness: 1.5pt, 
          paint: line_color.at(("a", "b", "c", "d", "e", "f", "g", "h").at(calc.rem(i, 8))),
          dash: (lines_style.b, lines_style.c, lines_style.d, lines_style.e).at(calc.rem(i, 4))
        ),
        mark: none,
        label: labels.at(i)
      )
    )
  }
  
  return plots
}
// 使用示例
#let filenames = ("M_R_B=20.csv", "M_R_B=15.csv", "M_R_B=10.csv", "M_R_B=5.csv")
#let labels = (
  [$B_"eff"=20 upright(M e V \/ f m^3)$],
  [$B_"eff"=15 upright(M e V \/ f m^3)$],
  [$B_"eff"=10 upright(M e V \/ f m^3)$],
  [$B_"eff"=5 upright(M e V \/ f m^3)$]
)

#let plots = generate_plots(filenames, labels)

// 然后将这些绘图命令作为参数传递给 diagram
#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $R (upright(k m))$  ],
    ylabel: [ $M \/ M_(dot.o)$],
    title: [ $"Einstein's Gravity", theta=0.0$ ],
    xlim: (0, 15.0),
    ylim: (-0.1, 2.5),
    xaxis:(subticks:1),
    yaxis:(subticks:3, ticks:range(0, 3, step:1)),
    legend: (position: top + left),
    ..plots  // 使用展开语法
  )
]

#let filenames = ("M_R_alpha=1.0e-6.csv", 
      "M_R_alpha=1.csv", 
      "M_R_alpha=3.csv", 
      "M_R_alpha=5.csv",
      "M_R_alpha=6.csv",
)
#let labels = (
  [$alpha=0.0$],
  [$alpha = 1 upright(k m)^2$],
  [$alpha = 3 upright(k m)^2$],
  [$alpha = 5 upright(k m)^2$],
  [$alpha = 6 upright(k m)^2$],
)

#let plots = generate_plots(filenames, labels)

// 然后将这些绘图命令作为参数传递给 diagram
#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $R (upright(k m))$  ],
    ylabel: [ $M \/ M_(dot.o)$],
    title: [ $"4-D EDB gravity", B_"eff"=10.0,theta=0.0$ ],
    xlim: (0, 15.0),
    ylim: (-0.1, 2.5),
    xaxis:(subticks:1),
    yaxis:(subticks:3, ticks:range(0, 3, step:1)),
    legend: (position: top + left),
    ..plots  // 使用展开语法
  )
]


#let filenames = (
  "M_R_theta=0.0.csv", 
  "M_R_theta=pi3.csv", 
  "M_R_theta=pi2.csv", 
  "M_R_theta=2pi3.csv",
  "M_R_theta=pi.csv",
)
#let labels = (
  [$theta=0.0$],
  [$theta=pi\/3$],
  [$theta=pi\/2$],
  [$theta=2pi\/3$],
  [$theta=pi$],
)

#let plots = generate_plots(filenames, labels)

// 然后将这些绘图命令作为参数传递给 diagram
#align(center)[
  #lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $R (upright(k m))$  ],
    ylabel: [ $M \/ M_(dot.o)$],
    title: [ $"4-D EDB gravity", B_"eff"=10.0,alpha = 6.0$ ],
    xlim: (0, 15.0),
    ylim: (-0.1, 2.5),
    xaxis:(subticks:1),
    yaxis:(subticks:3, ticks:range(0, 3, step:1)),
    legend: (position: top + left),
    ..plots  // 使用展开语法
  )
]