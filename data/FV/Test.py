from pypdf import PdfReader, PdfWriter

def split_pdf_pages(input_pdf, output_folder):
    """将PDF的每一页保存为单独的文件"""
    # 读取PDF文件
    reader = PdfReader(input_pdf)
    
    # 遍历每一页
    for page_num in range(len(reader.pages)):
        # 创建新的PDF写入器
        writer = PdfWriter()
        
        # 添加当前页
        writer.add_page(reader.pages[page_num])
        
        # 保存为单独的文件
        output_filename = f"{output_folder}/page_{page_num + 1}.pdf"
        with open(output_filename, 'wb') as output_file:
            writer.write(output_file)
    
    print(f"成功分割 {len(reader.pages)} 页")

# 使用示例
split_pdf_pages("plot.pdf", "output_pages")