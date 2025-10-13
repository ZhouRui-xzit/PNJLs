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


= Basic concepts of surface geometry

== First Fundamental Form and Second Fundamental Form

Let's consider  $2 d$ surface $S$ in $3 d$ Euclidean space. Thus,  $S$ can be 
parametrized by 
  $ 
    vb(r) = vb(r)(u,v) = (x(u,v), y(u,v), z(u,v)), (u,v) in U subset R^2. 
  $
For sphere：
  $ 
    x = R cos phi sin theta  \
    y = R sin phi sin theta  \
    z = R cos theta  \
    phi in [0,2pi), theta in [0,pi]  \
  $
whihc has been parametrized by spherical coordinates. Now,let   $ 
    &vb(r)_theta = pdv(vb(r), theta),
    vb(r)_phi = pdv(vb(r), phi)\  
    &vb(r)_(phi phi) = pdv(vb(r), phi,2),
    vb(r)_(theta theta) = pdv(vb(r), theta,2)\
    &vb(r)_(theta phi) =vb(r)_(phi theta)= pdv(vb(r), theta, phi)\
  $
Then，we have
  $ 
     vb(r)_theta &= (R cos phi cos theta, R sin phi cos theta, -R sin theta) \
      vb(r)_phi &= (-R sin phi sin theta, R cos phi sin theta, 0) \
  $
Now, we can get the first fundamental form of surface $S$:
  $ 
    E &= vb(r)_theta dot.c vb(r)_theta = R^2 \
    F &= vb(r)_theta dot.c vb(r)_phi = 0 \
    G &= vb(r)_phi dot.c vb(r)_phi = R^2 sin^2 theta \
  $
Now, we have the first fundamental form of surface $S$:
  $ 
    I = E dd(theta)^2 + 2F dd(theta)dd(phi) + G dd(phi)^2 =
    R^2 dd(theta^2) + R^2 sin^2 theta dd(phi)^2 \ 
  $
and the area element is 
  $ 
    dd(A)  = sqrt(E G-F^2) dd(theta)dd(phi) = R^2 sin theta dd(theta)dd(phi) \
  $
  Thus 
  $ 
    S = integral dd(A) = integral sqrt(E G-F^2) dd(theta)dd(phi)
    =
  $
  
 