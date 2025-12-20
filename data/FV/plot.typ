#import "@preview/lilaq:0.5.0" as lq
#import "@preview/codly:1.3.0": * // Code

#import "@preview/physica:0.9.5": *  // 物理符号和单位包
#import "@preview/numty:0.0.5" as nt  // 数值计算包

#import  "../plot_config.typ": *

// 
#set text(font: ("Libertinus Serif", "SimSun"), size: 12pt)


#show strong: text.with(font: ("Libertinus Serif", "SimHei"), size: 12pt)

#show emph: text.with(font: ("Libertinus Serif", "KaiTi"), size: 12pt)

#show math.equation: set text(font: ("Libertinus Math",), size: 12pt)




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


// 更好的方案：扩展 config 格式，添加 split 标记
#let load_cs_data(mu_B, R, configs) = {
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
    
    let base_path = "cs/rep_muB=" + repr(mu_B) + "_R=" + repr(R) + "_e=" + repr(e)
    let segments = ()
    
    if has_split {
      // 读取 _before 和 _after
      let before_data = csv(base_path + "_before.csv")
      let after_data = csv(base_path + "_after.csv")
      
      segments.push((
        T: before_data.slice(1).map(row => float(row.at(0))),
        P: before_data.slice(1).map(row => float(row.at(2))),
        E: before_data.slice(1).map(row => float(row.at(3))),
        TA: before_data.slice(1).map(row => float(row.at(4))),
        S: before_data.slice(1).map(row => float(row.at(5))),
        cs2: before_data.slice(1).map(row => float(row.at(7))),
        Cv: before_data.slice(1).map(row => float(row.at(8)))

      ))
      
      segments.push((
        T: after_data.slice(1).map(row => float(row.at(0))),
        P: after_data.slice(1).map(row => float(row.at(2))),
        E: after_data.slice(1).map(row => float(row.at(3))),
        TA: after_data.slice(1).map(row => float(row.at(4))),
        S: after_data.slice(1).map(row => float(row.at(5))),
        cs2: after_data.slice(1).map(row => float(row.at(7))),
        Cv: after_data.slice(1).map(row => float(row.at(8)))
      ))
    } else {
      // 读取 _all
      let data = csv(base_path + "_all.csv")
      
      segments.push((
        T: data.slice(1).map(row => float(row.at(0))),
        P: data.slice(1).map(row => float(row.at(2))),
        E: data.slice(1).map(row => float(row.at(3))),
        TA: data.slice(1).map(row => float(row.at(4))),
        S: data.slice(1).map(row => float(row.at(5))),
        cs2: data.slice(1).map(row => float(row.at(7))),
        Cv: data.slice(1).map(row => float(row.at(8)))
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
  (0.0, "b", "b", true),   // e=0.0, 虚线, 蓝色
  (0.3, "b", "b", true),   // e=0.3, 虚线, 红色
  (0.7, "c", "c"),   // e=0.7, 点线, 绿色
  (1.0, "d", "d"),   // e=1.0, 点划线, 橙色
)


#let cs_configs4 = (
  (0.0, "b", "a", true),   // e=0.0, 虚线, 蓝色
  (0.3, "b", "b", true),   // e=0.3, 虚线, 红色
  (0.7, "b", "c", true),   // e=0.7, 点线, 绿色
  (1.0, "d", "d", true),   // e=1.0, 点划线, 橙色
)



#let cs_data1 = load_cs_data(3.0, 30.0, cs_configs)
#let cs_data2 = load_cs_data(600.0, 30.0, cs_configs)
#let cs_data3 = load_cs_data(850.0, 30.0, cs_configs)
#let cs_data4 = load_cs_data(920.0, 30.0, cs_configs4)


//#let cs_data21 = load_cs_data(50.0, 100.0, cs_configs4)
//#let cs_data22 = load_cs_data(123.0, 100.0, cs_configs2)

//#let cs_data23 = load_cs_data(150.0, 100.0, cs_configs)

#let fit-to-width(width: 100%, content) = layout(size => context {
  let measured = measure(content)
  let target = if type(width) == ratio { size.width * width } else { width }
  let s = target / measured.width * 100%
  // 保持宽高比进行缩放
  scale(s, content)
})






#let fig_P1 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $P (upright(M e V \/ f m^3))$ ],
    title: [  $mu_B=0.0 upright(M e V)$  ],
    xlim: (0.0, 300.0),
    ylim: (0.0, 10000.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data1.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.P,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )



 #let fig_P2 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $P(upright(M e V \/ f m^3))$ ],
    title: [  $mu_B=600.0 upright(M e V)$  ],
    xlim: (10.0, 300.0),
    ylim: (0.0, 10000.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data2.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.P,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )





#let fig_P3 = lq.diagram(
    width: 10cm,
    height: 7cm,
     xlabel: [ $T$ (MeV) ],
    ylabel: [ $P (upright(M e V \/ f m^3))$ ],
    title: [  $mu_B=850.0 upright(M e V)$  ],
    xlim: (10.0, 300.0),
    ylim: (0.0, 10000.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data3.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.P,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )





#let fig_P4 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $P (upright(M e V \/ f m^3))$ ],
    title: [  $mu_B=920.0 upright(M e V)$  ],
    xlim: (10.0, 300.0),
    ylim: (0.0, 10000.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data4.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.P,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )




#let fig_E1 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $E (upright(M e V \/ f m^3))$ ],
    title: [  $mu_B=0.0 upright(M e V)$  ],
    xlim: (0.0, 300.0),
    ylim: (0.0, 100000.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data1.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.E,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )




#let fig_E2 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $E (upright(M e V \/ f m^3))$ ],
    title: [  $mu_B=600.0 upright(M e V)$  ],
    xlim: (0.0, 300.0),
    ylim: (0.0, 100000.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data2.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.E,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )




  #let fig_E3 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $E (upright(M e V \/ f m^3))$ ],
    title: [  $mu_B=850.0 upright(M e V)$  ],
    xlim: (0.0, 300.0),
    ylim: (0.0, 10000.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data3.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.E,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )



#let fig_E4 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $E (upright(M e V \/ f m^3))$ ],
    title: [  $mu_B=920.0 upright(M e V)$  ],
    xlim: (0.0, 300.0),
    ylim: (0.0, 10000.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data4.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.E,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )




#let fig_TA1 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $E-3P (upright(M e V \/ f m^3))$ ],
    title: [  $mu_B=0.0 upright(M e V)$  ],
    xlim: (0.0, 300.0),
    ylim: (0.0, 10000.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data1.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.TA,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )



#let fig_TA2 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $E-3P (upright(M e V \/ f m^3))$ ],
    title: [  $mu_B=600.0 upright(M e V)$  ],
    xlim: (0.0, 300.0),
    ylim: (0.0, 10000.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data2.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.TA,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )




#let fig_TA3 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $E-3P (upright(M e V \/ f m^3))$ ],
    title: [  $mu_B=850.0 upright(M e V)$  ],
    xlim: (0.0, 300.0),
    ylim: (0.0, 10000.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data3.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.TA,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )




#let fig_TA4 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $E-3P (upright(M e V \/ f m^3))$ ],
    title: [  $mu_B=920.0 upright(M e V)$  ],
    xlim: (0.0, 300.0),
    ylim: (0.0, 10000.0),
    yscale: "symlog",

    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data4.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.TA,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )



#let fig_S1 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $S (upright(f m^(-3)))$ ],
    title: [  $mu_B=0.0 upright(M e V)$  ],
    xlim: (0.0, 300.0),
    ylim: (0.0, 100.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data1.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.S,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )




#let fig_S2 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $S (upright(f m^(-3)))$ ],
    title: [  $mu_B=600.0 upright(M e V)$  ],
    xlim: (0.0, 300.0),
    ylim: (0.0, 100.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data2.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.S,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )





#let fig_S3 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $S (upright(f m^(-3)))$ ],
    title: [  $mu_B=850.0 upright(M e V)$  ],
    xlim: (0.0, 300.0),
    ylim: (0.0, 100.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data3.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.S,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )



#let fig_S4 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $S (upright(f m^(-3)))$ ],
    title: [  $mu_B=920.0 upright(M e V)$  ],
    xlim: (0.0, 300.0),
    ylim: (0.0, 100.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data4.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.S,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )

#let fig_CV1 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $C_V$ ],
    title: [  $mu_B=0.0 upright(M e V)$  ],
    xlim: (0.0, 300.0),
    ylim: (0.0, 500.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data1.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.Cv,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )


#let fig_CV2 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $C_V$ ],
    title: [  $mu_B=600.0 upright(M e V)$  ],
    xlim: (0.0, 300.0),
    ylim: (0.0, 500.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data2.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.Cv,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )


#let fig_CV3 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $C_V$ ],
    title: [  $mu_B=850.0 upright(M e V)$  ],
    xlim: (0.0, 300.0),
    ylim: (0.0, 500.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data3.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.Cv,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )




#let fig_CV4 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $C_V$ ],
    title: [  $mu_B=920.0 upright(M e V)$  ],
    xlim: (0.0, 300.0),
    ylim: (-0.0, 500.0),
    yscale: "symlog",
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data4.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.Cv,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )






#let fig_cs1 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $c_(s\/rho_B)$ ],
    title: [  $mu_B=0.0 upright(M e V)$  ],
    xlim: (50.0, 300.0),
    ylim: (-0.0, 0.35),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data1.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.cs2,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )




#let fig_cs2 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $c_(s\/rho_B)$ ],
    title: [  $mu_B=600.0 upright(M e V)$  ],
    xlim: (50.0, 300.0),
    ylim: (-0.0, 0.35),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data2.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.cs2,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )




#let fig_cs3 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $c_(s\/rho_B)$ ],
    title: [  $mu_B=850.0 upright(M e V)$  ],
    xlim: (50.0, 300.0),
    ylim: (-0.0, 0.35),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data3.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.cs2,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )



#let fig_cs4 = lq.diagram(
    width: 10cm,
    height: 7cm,
    xlabel: [ $T$ (MeV) ],
    ylabel: [ $c_(s\/rho_B)$ ],
    title: [  $mu_B=920.0 upright(M e V)$  ],
    xlim: (50.0, 300.0),
    ylim: (-0.0, 0.35),
    xaxis:(subticks:1),
    yaxis:(subticks:1),
    legend: (position: top + left),
    
    ..cs_data4.values().map(d => {
      d.segments.enumerate().map(((i, seg)) => 
        lq.plot(
          seg.T, seg.cs2,
          stroke: (thickness: 1.2pt, paint: d.color, dash: d.style, cap: "round"),
          mark: none,
          label: if i == 0 { [ $delta=#d.e$ ] } else { none }
        )
      )
    }).flatten()
  )



#pagebreak()
// 2. 在表格中使用
#table(
  columns: (auto, auto),
  inset: 10pt,
  align: horizon+center,
  stroke: 1pt + gray,
  [#fig_P1],[#fig_P2],
  [#fig_P3],[#fig_P4],
 
)
#pagebreak()

#table(
  columns: (auto, auto),
  inset: 10pt,
  align: horizon+center,
  stroke: 1pt + gray,
  [#fig_E1],[#fig_E2],
  [#fig_E3],[#fig_E4],
 
)
#pagebreak()
#table(
  columns: (auto, auto),
  inset: 10pt,
  align: horizon+center,
  stroke: 1pt + gray,
  [#fig_TA1],[#fig_TA2],
  [#fig_TA3],[#fig_TA4],
 
)

#pagebreak()
#table(
  columns: (auto, auto),
  inset: 10pt,
  align: horizon+center,
  stroke: 1pt + gray,
  [#fig_S1],[#fig_S2],
  [#fig_S3],[#fig_S4],
 
)
#pagebreak()

#table(
  columns: (auto, auto),
  inset: 10pt,
  align: horizon+center,
  stroke: 1pt + gray,
  [#fig_CV1],[#fig_CV2],
  [#fig_CV3],[#fig_CV4],
 
)

#pagebreak()

#table(
  columns: (auto, auto),
  inset: 10pt,
  align: horizon+center,
  stroke: 1pt + gray,
  [#fig_cs1],[#fig_cs2],
  [#fig_cs3],[#fig_cs4],
 
)