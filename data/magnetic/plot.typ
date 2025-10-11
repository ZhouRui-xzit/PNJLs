#import "@preview/lilaq:0.5.0" as lq
#import "@preview/codly:1.3.0": * // Code

#import "@preview/physica:0.9.5": *  // 物理符号和单位包
#import "@preview/numty:0.0.5" as nt  // 数值计算包


// 设置正文字体为Arial（英文）和方正书宋（中文）
#set text(font: ("Arial", "FZShuSong-Z01"), size: 12pt)

// 设置粗体字体为Arial（英文）和方正黑体（中文）
#show strong: text.with(font: ("Arial","FZHei-B01"), size: 12pt)

// 设置斜体字体为Arial（英文）和方正楷体（中文）
#show emph: text.with(font: ("Arial", "FZKai-Z03"), size: 12pt)

// 设置数学公式字体
#show math.equation: set text(font: "Latin Modern Math")








#let line_color = (
  a: rgb("#0b86eb"),    // 深蓝
  b: rgb("#ea074fb6"),     // 砖红
  c: rgb(44, 160, 44),     // 翠绿
  d: oklab(70.63%, 0.111, 0.141),    // 橙色
  e: rgb("#7e0aea")    // 淡紫
)

#let lines_style = (
  a: "solid",                      // 实线
  b: "dashed",                     // 虚线
  c: "dotted",                     // 点线
  d: "dash-dotted",                // 点划线
  e: (5pt, 2pt),                   // 长虚线
  f: (2pt, 2pt),                   // 短虚线
  g: (3pt, 0.5pt, 1pt, 0.5pt),         // 密集点划
  h: (6pt, 3pt, 1pt, 3pt),         // 密集点线
)

// 使用示例:
// stroke: (thickness: 1pt, paint: colors.a, dash: lines.b)




= mu-rho 曲线

#let SSS = lq.load-txt(read("T=20.csv"), skip-rows:1)

#let S1 = SSS.at(0)
#let S2 = SSS.at(1)
#let S3 = SSS.at(2)

#align(center)[
  #lq.diagram(
  width: 10cm,
  height: 7cm,
  xlabel: [ $mu_B$ (MeV) ],
  ylabel: [ $rho_B\/rho_0$ ],
  title: [ $mu_B-rho_q$ ],
  xlim: (0.0, 3.0),
  ylim: (200, 400),
  xaxis:(subticks:1),
  yaxis:(subticks:1, ticks:range(200, 401, step:50)),
  lq.plot(
    S2, S3,
    stroke: (thickness:1pt, paint: blue, dash: lines_style.d),
    mark: none,  // 不显示散点标记
    //smooth: true,

  ),
)
]

= phi-T 曲线

#let T1 = lq.load-txt(read("Tmu_eB01.csv"), skip-rows:1)
#let T2 = lq.load-txt(read("Tmu_eB02.csv"), skip-rows:1)
#let T3 = lq.load-txt(read("Tmu_eB03.csv"), skip-rows:1)
#let T4 = lq.load-txt(read("Tmu_eB05.csv"), skip-rows:1)

#let Ts = T1.at(0)
#let phiu1 = T1.at(1)
#let phiu2 = T2.at(1)
#let phiu3 = T3.at(1)
#let phiu4 = T4.at(1)

#align(center)[
  #lq.diagram(
  width: 10cm,
  height: 7cm,
  xlabel: [ $T$ (MeV) ],
  ylabel: [ $phi.alt_u$ ],
  title: [ $T-phi.alt_u$ ],
  xlim: (50.0, 300.0),
  ylim: (-3.0, 0.0),
  xaxis:(subticks:1),
  yaxis:(subticks:1, ticks:range(-3,1, step:1)),
  legend: (position: top + left),
  lq.plot(
    Ts, phiu1,
    stroke: (thickness:1.5pt, paint: line_color.a, dash:lines_style.g),
    mark: none,
    label:[ $e B = 0.1 upright(G e V)$ ]
  ),
  lq.plot(
    Ts, phiu2,
    stroke: (thickness:1.5pt, paint: line_color.b, dash:lines_style.h),
    mark: none,
    label:[ $e B = 0.2 upright(G e V)$ ]
  ),
  lq.plot(
    Ts, phiu3,
    stroke: (thickness:1.5pt, paint: line_color.c, dash:lines_style.c),
    mark: none,
    label:[ $e B = 0.3 upright(G e V)$ ]
  ),
  lq.plot(
    Ts, phiu4,
    stroke: (thickness:1.5pt, paint: line_color.d, dash:lines_style.b),
    mark: none,
    label:[ $e B = 0.5 upright(G e V)$ ]
),
  )]



#let T_mu_all_eB = csv("Tmu_all_eB.csv")  // 读取CSV文件，包含表头
#let T_mu_all_eB_num = T_mu_all_eB.slice(1).map(row => row.map(x => float(x))) 
// 转换为全float的二维数组（不含header）

#let filter_by_eB(arr, value) = arr.filter(row => row.at(1) == value) // 函数：根据eB值筛选数据行



 // 计算并返回指定eB值的数据
#let calc_sigma_data(eB_value) = {
  // 筛选指定eB值的所有数据行
  let rows = filter_by_eB(T_mu_all_eB_num, eB_value)
  
  // 提取T, sigma_u和sigma_d列
  let Ts = rows.map(row => row.at(0))
  let sigma_u = rows.map(row => row.at(2))
  let sigma_d = rows.map(row => row.at(3))
  
  // 计算sigma_plus和sigma_minus
  let sigma_plus = nt.mult(nt.add(sigma_u, sigma_d), 0.5)
  let sigma_minus = nt.sub(sigma_u, sigma_d)
  
  // 返回结果字典
  (
    T: Ts,
    sigma_u: sigma_u,
    sigma_d: sigma_d,
    sigma_plus: sigma_plus,
    sigma_minus: sigma_minus,
    eB: eB_value
  )
}

// 使用示例
#let Sigma_eB02 = calc_sigma_data(0.2)
#let Sigma_eB04 = calc_sigma_data(0.4)
#let Sigma_eB06 = calc_sigma_data(0.6)
#let Sigma_eB08 = calc_sigma_data(0.8)


#let fig1(w, h) = align(center)[
  #lq.diagram(
  width: w,
  height: h,
  xlabel: [ $T$ (MeV) ],
  ylabel: [ $(Sigma_u+Sigma_d)\/2$ ],
  //title: [ $T-phi.alt_u$ ],
  xlim: (50.0, 300.0),
  ylim: (-0.1, 2.00),
  xaxis:(subticks:1),
  yaxis:(subticks:1, ticks:range(-1, 3, step:1)),
  legend: (position: top + right),
  lq.plot(
    Sigma_eB02.T, Sigma_eB02.sigma_plus,
    stroke: (thickness:1pt, paint: blue),
    mark: none,
    label:[ $e B = 0.2 upright(G e V)$ ]
  ),
  lq.plot(
    Sigma_eB04.T, Sigma_eB04.sigma_plus,
    stroke: (thickness:1pt, paint: red),
    mark: none,
    label:[ $e B = 0.4 upright(G e V)$ ]
  ),
  lq.plot(
    Sigma_eB06.T, Sigma_eB06.sigma_plus,
    stroke: (thickness:1pt, paint: green),
    mark: none,
    label:[ $e B = 0.6 upright(G e V)$ ]
  ),
  lq.plot(
    Sigma_eB08.T, Sigma_eB08.sigma_plus,
    stroke: (thickness:1pt, paint: orange),
    mark: none,
    label:[ $e B = 0.8 upright(G e V)$ ]
  ),
  )]


#let fig2(w,h)=align(center)[
  #lq.diagram(
  width: w,
  height: h,
  xlabel: [ $T$ (MeV) ],
  ylabel: [ $(Sigma_u-Sigma_d)$ ],
  //title: [ $T-phi.alt_u$ ],
  xlim: (50.0, 300.0),
  ylim: (-0.1, 1.00),
  xaxis:(subticks:1),
  yaxis:(subticks:1, ticks:range(-1, 3, step:1)),
  legend: (position: top + right),
  lq.plot(
    Sigma_eB02.T, Sigma_eB02.sigma_minus,
    stroke: (thickness:1pt, paint: blue),
    mark: none,
    label:[ $e B = 0.2 upright(G e V)$ ]
  ),
  lq.plot(
    Sigma_eB04.T, Sigma_eB04.sigma_minus,
    stroke: (thickness:1pt, paint: red),
    mark: none,
    label:[ $e B = 0.4 upright(G e V)$ ]
  ),
  lq.plot(
    Sigma_eB06.T, Sigma_eB06.sigma_minus,
    stroke: (thickness:1pt, paint: green),
    mark: none,
    label:[ $e B = 0.6 upright(G e V)$ ]
  ),
  lq.plot(
    Sigma_eB08.T, Sigma_eB08.sigma_minus,
    stroke: (thickness:1pt, paint: orange),
    mark: none,
    label:[ $e B = 0.8 upright(G e V)$ ]
  ),
  )]


#align(center)[
  #table(
  columns: (auto, auto),

  inset: 10pt,
  align: horizon+center,
  stroke:  1pt + gray,
  table.header(
    [ $ (Sigma_u+Sigma_d)\/2 $ ],[ $ (Sigma_u-Sigma_d) $ ]
  ),
  box(inset: 1pt)[#scale(65%, fig1(10cm, 6cm))], 
    box(inset: 1pt)[#scale(65%, fig2(10cm, 6cm))]
  )
]