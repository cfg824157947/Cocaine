import gseapy as gp
import pandas as pd

# Load your data (example: gene expression dataset and rank file)
# Ensure the data has two columns: 'gene' and 'score'
module_dir="/project/zduren/durenlab/palmetto/cham/Heroin/LINGER/Module_genes/"
#module_file = module_dir + "K15p0.7000000000000001M1.csv"
p_list=[0.5]
#p_list=[0.5,0.55,0.6]
for p in p_list:
    module_file = module_dir + "K10p"+str(p) + "M1.csv"


    module_genes = pd.read_csv(module_file, index_col=0)
    #module_genes=pd.read_csv()
    DEG_dir = "/project/zduren/durenlab/palmetto/cham/Heroin/analyze/NInteResult/RNA/PROTEIN/3/Gene_score/"
    entry = DEG_dir + "L6_IT_Car3.csv"
    #entry = DEG_dir + "L5_IT.csv"
    #entry = DEG_dir + "OPC.csv"
    for entry in os.listdir(DEG_dir):
        celltype = entry.split('.')[0]
        DEG_file = os.path.join(DEG_dir, entry)
        ranked_genes = pd.read_csv(DEG_file, index_col=2)  # Replace with your ranked list file
        ranked_genes = ranked_genes[['t']]  # Replace with your ranked list file
        #ranked_genes = ranked_genes[['logFC']]  # Replace with your ranked list file
        ranked_genes = ranked_genes.loc[ranked_genes.index.isin(module_genes.index)]
        # Example format:
        # gene,score
        # GeneA,2.5
        # GeneB,-1.3

        # Perform GSEA analysis
        gsea_result_dir = "/project/zduren/durenlab/palmetto/cham/Heroin/LINGER/gsea_results/"
        output_dir = os.path.join(gsea_result_dir + celltype + '/')
        if not os.path.exists(output_dir):
            # Create the directory, including intermediate directories
            os.makedirs(output_dir, exist_ok=True)
            print(f"Directory '{output_dir}' was created.")

        gsea_results = gp.prerank(
            rnk=ranked_genes,  # Ranked gene list
            gene_sets='KEGG_2021_Human',  # Predefined gene sets; you can also use custom gene sets
            outdir=output_dir,  # Directory for results
            permutation_num=1000,  # Number of permutations for p-value calculation
        #    min_size=15,  # Minimum size of gene sets
        #    max_size=500,  # Maximum size of gene sets
            seed=123,  # Reproducibility seed
        )

        # View and save results
        #gsea_results.res2d.to_csv("gsea_report.csv")
#        print(gsea_results.res2d)  # Shows the results summary

        figure_dir =  "/project/zduren/durenlab/palmetto/cham/Heroin/LINGER/figure/GSEA/"
        figure_dir =  "/project/zduren/durenlab/palmetto/cham/Heroin/LINGER/figure/GSEA/"
        term_list = ["Morphine addiction"]
        for term in term_list:
            figure_file = 'figure/GSEA/'+ celltype+ str(p) +term+'.png'
       # gsea_results.plot(
       #     terms="MAPK signaling pathway",  # Replace with the gene set name you want to plot
       #     ofname='figure/GSEA/MAPK_gsea_results.png'# Save as a PNG file
       # )
            gsea_results.plot(
                terms=term,  # Replace with the gene set name you want to plot
                ofname=figure_file# Save as a PNG file
            )


