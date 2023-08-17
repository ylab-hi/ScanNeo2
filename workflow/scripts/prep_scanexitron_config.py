import sys
import os

def main():
    template = """
    [fasta]

    hg38={hg38_path}
    hg19=

    [annotation]

    hg38={hg38_annotation_path}
    hg19=

    [cds]
    
    hg38={hg38_cds_path}
    hg19=
    """
    parameters = {
            'hg38_path': str(os.path.abspath(sys.argv[1])),
            'hg38_annotation_path': str(os.path.abspath(sys.argv[2])),
            'hg38_cds_path': str(os.path.abspath(sys.argv[3]))
    }

    filled_template = template.format(**parameters)
    output_path = sys.argv[4]

    with open(output_path, 'w') as output_file:
        output_file.write(filled_template)

main()
