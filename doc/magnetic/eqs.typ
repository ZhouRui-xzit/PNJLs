//==========================================
// 1. 基础包导入
//==========================================
#import "@preview/physica:0.9.5": *  // 物理符号和单位包
#import "@preview/hydra:0.6.2": *    // 多功能样式包
#import "@preview/codly-languages:0.1.1": *
#import "@preview/ctheorems:1.1.3": * //Box

#import "@preview/cetz:0.4.2"        // 示意图功能包
#import "@preview/lilaq:0.4.0" as lq // 科学绘图功能包
#import "@preview/codly:1.3.0": * // Code


//==========================================
// 2. 文档基本设置
//==========================================
#set document(title: "Theoretical Physics Research", author: "Author Name")

//==========================================
// 3. 页面布局设置
//==========================================
#set page(
  paper: "a4",
  fill: white, 
  margin: (x:40pt, y:60pt), 
  numbering: "1", 
  // 页眉设置：右对齐蓝色hydra图标，下方添加分隔线
  header: context {
    align(right, text(size: 13pt, fill: blue, hydra(2)))   
    v(-10pt)
    line(length: 100%, stroke: black) 
  }, 
  // 页脚设置：居中页码，上方添加分隔线
  footer: [
    #align(center, context {counter(page).display("1")})
    #line(length: 100%, stroke: black) 
  ]
)

// 设置脚注文本为红色
#show footnote.entry: set text(fill: rgb("#f11809"))

// 添加目录并设置标题
#outline(indent: auto, title: "Table of Contents")
#pagebreak()

// 重置页码从1开始（目录后的第一页）
#counter(page).update(1)

//==========================================
// 4. 标题样式设置
//==========================================
// 设置标题编号格式为"1.1."
#set heading(numbering: "1.1.")

// 为所有标题添加下方空间
#show heading: it => {
  it
  par()[#text(size:0.5em)[#h(0.0em)]]
}

// 为所有图表添加下方空间
#show figure: it => {
  it
  par()[#text(size:0.5em)[#h(0.0em)]]
}

// 一级标题：添加分页、居中、蓝色、加粗、18pt
#show heading.where(level: 1): it => pagebreak(weak: true) + it
#show heading.where(level: 1): it => {
  align(center, text(
    size: 18pt,
    fill: blue,
    weight: "bold",
    it
  ))
  v(-1pt)
}

// 二级标题：蓝色、16pt
#show heading.where(level: 2): it => {
  text(
    size: 16pt,
    fill: blue,
    it
  )
  v(-10pt)
}

// 三级标题：蓝色、14pt
#show heading.where(level: 3): it => {
  text(
    size: 14pt,
    fill: blue,
    it
  )
  v(-10pt)
}

// 设置默认对齐方式为左对齐+顶部对齐
#set align(left+top)

//==========================================
// 5. 文本字体设置
//==========================================
// 设置正文字体为Arial（英文）和方正书宋（中文）
#set text(font: ("Arial", "FZShuSong-Z01"), size: 12pt)

// 设置粗体字体为Arial（英文）和方正黑体（中文）
#show strong: text.with(font: ("Arial","FZHei-B01"), size: 12pt)

// 设置斜体字体为Arial（英文）和方正楷体（中文）
#show emph: text.with(font: ("Arial", "FZKai-Z03"), size: 12pt)

// 设置数学公式字体
#show math.equation: set text(font: "Latin Modern Math")

//==========================================
// 6. 图表和公式编号设置
//==========================================
// 设置图表编号格式：(章节号.序号)
#set figure(numbering: (..nums) => {
  let ch = counter(heading).get().first()
  numbering("(1.1)", ch, ..nums)
})



// 设置公式编号格式：(两层)

#set math.equation(numbering: (..nums) => {
  let ch = counter(heading).get().first()
  numbering("(1.1)", ch, ..nums)
})


// 设置公式、表格和图像的前缀
#set math.equation(supplement: [式])

// 设置表格前缀和标题位置
#show figure.where(kind: table): set figure(supplement: [表])
#show figure.where(kind: table): set figure.caption(position: top)

// 设置图像前缀
#show figure.where(kind: image): set figure(supplement: [图])

// 二级标题后重置公式编号
#show heading.where(level: 1): it => it + counter(math.equation).update(0)

// 一级标题后重置图像和表格编号
#show heading.where(level: 1): it => it + counter(figure.where(kind:image)).update(0)
#show heading.where(level: 1): it => it + counter(figure.where(kind:table)).update(0)

//==========================================
// 7. 引用和链接样式设置
//==========================================
// 设置链接文本为红色
#show link: text.with(fill: red)

// 设置引用文本为红色
#show ref: it => {
  text(it, red, font: ("Arial", "FZShuSong-Z01"), size: 12pt)
}

// 处理无效引用，显示红色错误提示
#show ref: it => {
  if query(it.target).len() == 0 {
    return text(fill: red, "<???" + ">")
  }
  it
}



#let thm = thmbox("theorem", "Theorem", 
fill: rgb("#e6d2b8"),inset: (x: 1.2em, top: 1em, bottom: 1em),
base_level: 1, supplement:"定理").with(numbering: "1.1")

#let def = thmbox("def", "Definition", 
fill:cmyk(30.61%, 1.22%, 0%, 3.92%),
inset: (x: 1.2em, top: 1em, bottom: 1em),base_level: 1,
supplement:"定义").with(numbering: "1.1")

#let exm = thmbox("exm", "Example", 
fill:rgb("#afdbb8"),
inset: (x: 1.2em, top: 1em, bottom: 1em),base_level: 1,
supplement:"例子").with(numbering: "1.1")

#let proposition = thmbox("Proposition", "Proposition", 
fill:color.hsv(196.19deg, 84.75%, 87.45%, 56.8%),inset: (x: 1.2em, top: 1em, bottom: 1em),base_level: 1,supplement:"命题").with(numbering: "1.1")

#let axm = thmbox("axiom", "Axiom", fill:rgb("#f3dfc5"),inset: (x: 1.2em, top: 1em, bottom: 1em),base_level: 1,supplement:"公理").with(numbering: "1.1")



#let proof = thmplain(
"proof",
"Proof",
bodyfmt: body => [
#text(emph(body)) #h(1fr) $qed$ 
],
titlefmt: it => {
text(emph(it))
}
).with(numbering: none)


#let remark = thmplain(
  "remark",
  "Remark",
  bodyfmt: it => [
    #text(emph(it)) 
    #h(1fr)
    $qed$ 
  ],
  titlefmt: it => {text(emph(it))},
  ).with(numbering: none)


= p-NJL in magnetic field

  $ 
    Omega(T,mu_B) = 2 G sum_f phi.alt_f^2 - 4 K phi.alt_u phi.alt_d phi.alt_s 
    + U(Phi,overline(Phi),T) + sum_f 
    (Omega_f^0 + Omega^"mag"_f + Omega_f^"T")  
  $
where  $ 
     Omega_f^0 &= -6 integral^Lambda dd(p,3)/(2pi)^3 
     sqrt(vb(p)^2+M_f^2) \ 
      Omega_f^"mag" &=  (-3 (abs(q_f) e B)^2)/(2pi^2) 
    [zeta'(-1,x_f)-1/2 (x_f^2-x_f) ln x_f +x_f^2/4] \ 
    Omega_f^"T" &=-T (abs(q_f) e B)/(2pi) 
    sum_(n=0)^infinity dd(p_z)/(2pi) 
    (Z_f^+ (E_f) + Z_f^- (E_f)) 
  $
and  $ 
    M_f &= m_(0 f) - 4 G phi.alt_f +2 K phi.alt_f' phi.alt_f'' \ 
  Z_f^+ (E_f) &= ln (1+3 Phi e^(-beta (E_f-mu))+ 3 overline(Phi) e^(-2beta (E_f-mu))
  + e^(-3 beta (E_f-mu))
  ) \ 
  Z_f^+ (E_f) &= ln (1+3overline(Phi)   e^(-beta (E_f+mu))+ 3 Phi e^(-2beta (E_f+mu))
  + e^(-3 beta (E_f+mu))
  ) \ 
  E_f &= sqrt(p_z^2+M_f^2+ 2n abs(q_f) e B)\ 
  x_f &= M_f^2/(2 e B abs(q_f))
  $
In particular:
  $ 
     zeta'(-1,x) = evaluated(dv(, s) zeta(s,x))_(s=-1)
  $
and 

  
For vacumun part, we have:
  $ 
    integral_0^Lambda p^2 sqrt(p^2+M_f^2) dd(p) 
    = 1/8 [Lambda sqrt(Lambda^2 +M_f^2) (2Lambda^2+M_f^2)
    - M_f^4 ln ((Lambda+sqrt(Lambda^2+M_f^2))/M_f)
    ]  
  $
  
i.e.  $ 
    Omega_f^0 = -3/(8 pi^2)   [Lambda sqrt(Lambda^2 +M_f^2) (2Lambda^2+M_f^2)
    - M_f^4 ln ((Lambda+sqrt(Lambda^2+M_f^2))/M_f)
    ]  
  $
  
  
  
  
  
  