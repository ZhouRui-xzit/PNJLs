#import "@preview/lilaq:0.5.0" as lq
#import "@preview/codly:1.3.0": * // Code

#import "@preview/physica:0.9.5": *  // 物理符号和单位包
#import "@preview/numty:0.0.5" as nt  // 数值计算包



#let line_color = (
  a: rgb("#0b86eb"),    // 深蓝
  b: rgb("#ea074fb6"),     // 砖红
  c: rgb(44, 160, 44),     // 翠绿
  d: oklab(70.63%, 0.111, 0.141),    // 橙色
  e: rgb("#7e0aea"),    // 淡紫
  f: rgb("#e377c2"), // 粉红色
  g: rgb("#8c564b"), // 棕色
  h: rgb("#17becf"), // 青色
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


// 变参数绘制phi T 图像 ,参数第二列，T 第一列 phi第三列.....


#let get_all_R_values(path) = {
  // 读取CSV文件并转换为浮点数数组
  let data = csv(path).slice(1).map(row => row.map(x => float(x)))
  
  // 获取所有R值（第二列）
  let all_R = data.map(row => row.at(1))
  
  // 创建唯一值集合
  let unique_R = ()
  for r in all_R {
    // 检查是否已存在于结果数组中
    let exists = false
    for ur in unique_R {
      if calc.abs(r - ur) < 1e-6 {
        exists = true
        break
      }
    }
    
    // 如果不存在，则添加
    if not exists {
      unique_R.push(r)
    }
  }
  

  return unique_R
}



#let fig_phiT(path) = {
  
  let values = get_all_R_values(path)
  let data = csv(path).slice(1).map(row => row.map(x => float(x)))
  
  // 创建结果字典
  let result = (:)
  

  // 为每个value值收集对应的数据
  for value in values {
    // 使用浮点数精度容差进行过滤
    let epsilon = 1e-6
    let filtered_rows = data.filter(row => {
      let diff = calc.abs(row.at(1) - value)
      diff < epsilon
    })
    
    if filtered_rows.len() > 0 {
      // 提取T和phi列
      let Ts = filtered_rows.map(row => row.at(0))
      let phis = filtered_rows.map(row => row.at(2))
      
      // 将这对值添加到结果中，使用value作为键
      result.insert(str(value), (T: Ts, phi: phis))
    } else {
      // 如果找不到值，添加空数据
      result.insert(str(value), (T: (), phi: ()))
    }
  }
  
  return result
}
