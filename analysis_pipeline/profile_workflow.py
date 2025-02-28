import graphviz
import os

def make_png(input_file_name, output_file_name):
    dot = graphviz.Source.from_file(input_file_name)
    dot.render(outfile=output_file_name)

if __name__ == '__main__':
    os.chdir(r"C:\Program Files\Graphviz\bin")
    make_png(r"C:\Users\spejo\Documents\1_CRISPR_analysis_test_input\FASTQ\download_fastq_files_baseline.dot", r"C:\Users\spejo\Documents\1_CRISPR_analysis_test_input\FASTQ\download_fastq_files_baseline.jpg")