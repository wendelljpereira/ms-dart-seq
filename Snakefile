configfile: "config.yaml"
run_in_slurm_env = True

if run_in_slurm_env == True:
    shell.prefix("module load R/3.5.1; module load trimmomatic/0.36; module load bowtie2/2.3.5; module load samtools/1.8; module load bedtools/2.27.1; module load subread/1.6.2; module load fastqc/0.11.4; module load alfa/1.1.1; ")

rule all:
    input:

rule barcode_removing:
    input:
        expand("{barcodes_files}", barcodes_files=config["barcodes_files"]),
        expand("{adapters_file}", adapters_file=config["adapters_file"]),
        expand("{fastq}.FASTQ", fastq=config["fastq"])
    output:
        expand("{fastq}.without_barcodes_and_adapters.fq", fastq=config["fastq"])
    shell:
        """
        # defines the input files

        barcode_file1={input[0]}
        barcode_file2={input[1]}
        adapters_file={input[2]}

        # execute using the two sets of barcodes

        while read FILE barcode_sequence barcode_size
        do
        trimmomatic SE -threads 12 $FILE.FASTQ $FILE.without_barcodes_and_adapters.fq ILLUMINACLIP:$adapters_file:2:30:6 HEADCROP:$barcode_size MINLEN:20
        done < "$barcode_file1"

        while read FILE barcode_sequence barcode_size
        do
        trimmomatic SE -threads 12 $FILE.FASTQ $FILE.without_barcodes_and_adapters.fq ILLUMINACLIP:$adapters_file:2:30:6 HEADCROP:$barcode_size MINLEN:20
        done < "$barcode_file2"

        """

rule fastqc:
    input:
        expand("{fastq}.without_barcodes_and_adapters.fq", fastq=config["fastq"])
    output:
        expand("fastqc_after_trimming/{fastq}.without_barcodes_and_adapters_fastqc.html", fastq=config["fastq"]),
        expand("fastqc_after_trimming/{fastq}.without_barcodes_and_adapters_fastqc.zip", fastq=config["fastq"])
    shell:
        """
        fastqc *.fq
        mv *fastqc.zip fastqc_after_trimming
        mv *fastqc.html fastqc_after_trimming
        """

rule combine_samples:
    input:
        expand("{fastq}.without_barcodes_and_adapters.fq", fastq=config["fastq"])
    output:
        expand("{combine_samples}", combine_samples=config["combine_samples"])
    shell:
        "cat {input} > {output}"

rule mapping:
    input:
        expand("{reference_genome}", reference_genome=config["reference_genome"]),
        expand("{combine_samples}", combine_samples=config["combine_samples"])
    output:
        expand("mapping/{mapping}.fq", mapping=config["mapping"]["fastq"]),
        expand("mapping/{mapping}.sam", mapping=config["mapping"]["sam"])
    shell:
        """

            base=$(basename {input[0]} ".fa")

            bowtie2-build -f --threads 24 {input[0]} $base

            bowtie2 -p 24 -q -k 2 -N 0 -L 15 -i S,1,0.50 -R 3 --no-1mm-upfront --end-to-end --no-unal --score-min L,0,0 -x $base --un {output[0]} -U {input[1]} -S {output[4]}

            bowtie2 -p 24 -q -k 2 -N 0 -L 15 -i S,1,0.50 -R 3 --no-1mm-upfront --end-to-end --no-unal --mp 2 --np 10 --rdg 10 --rfg 10 --score-min L,-2,0 -x $base --un {output[1]} -U {output[0]} -S {output[5]}

            bowtie2 -p 24 -q -k 2 -N 0 -L 15 -i S,1,0.50 -R 3 --no-1mm-upfront --end-to-end --no-unal --mp 2 --np 10 --rdg 10 --rfg 10 --score-min L,-4,0 -x $base --un {output[2]} -U {output[1]} -S {output[6]}

            bowtie2 -p 24 -q -k 2 -N 0 -L 15 -i S,1,0.50 -R 3 --no-1mm-upfront --end-to-end --no-unal --mp 2 --np 10 --rdg 10 --rfg 10 --score-min L,-6,0 -x $base --un {output[3]} -U {output[2]} -S {output[7]}

        """

rule extract_sing_mapping:
    input:
        expand("mapping/{mapping}.sam", mapping=config["mapping"]["sam"])
    output:
        expand("mapping/{mapping}.single_map.sam", mapping=config["mapping"]["sam"])
    shell:
        """

            for sample in `ls mapping/*.sam`
            do
            base=$(basename $sample \".sam\")

            samtools view -h mapping/${{base}}.sam | grep -v \"XS:i:\" > mapping/${{base}}.single_map.sam
            done

        """

rule sam_to_bam:
    input:
        expand("mapping/{mapping}.single_map.sam", mapping=config["mapping"]["sam"]),
    output:
        expand("mapping/{mapping}.single_map.bam", mapping=config["mapping"]["sam"])
    shell:
        """

            for sample in `ls mapping/*single_map.sam`
            do
            base=$(basename $sample \".sam\")
            samtools view -b mapping/${{base}}.sam > mapping/${{base}}.bam
            samtools sort mapping/${{base}}.bam > mapping/tmp.bam
            mv mapping/tmp.bam mapping/${{base}}.bam
            done

        """

rule merge_bam:
    input:
        expand("mapping/{mapping}.single_map.bam", mapping=config["mapping"]["sam"])
    output:
        expand("mapping/{combine_bam}", combine_bam=config["combine_bam"])
    shell:
        """

            samtools merge -f {output[0]} {input[0]} {input[1]} {input[2]} {input[3]}
            samtools sort {output[0]} > mapping/temp.bam

            mv mapping/temp.bam {output[0]}

        """

rule bam_to_bed:
    input:
        expand("mapping/{combine_bam}", combine_bam=config["combine_bam"])
    output:
        expand("bed_files/{bed_name}", bed_name=config["bam_to_bed"]["bed_name"])
    shell:
        """
            #conversion from bam to bed

            #Deletes intermediate files

            #rm mapping/*.sam
            #rm mapping/*map.bam

            bamToBed -i {input[0]} > {output[0]}
            bedtools sort -i {output[0]} > tmp.bed
            mv tmp.bed {output[0]}

        """

rule sample_site_definition:
    input:
        expand("bed_files/{bed_name}", bed_name=config["bam_to_bed"]["bed_name"]),
        expand("{enzymes_sites}", enzymes_sites=config["enzymes_sites"])
    output:
        expand("sample_sites_generation/{sample_site_definition}", sample_site_definition=config["sample_site_definition"])
    shell:
        """

            cat {input[0]} {input[1]} > bed_files/enz_sites_sample_reads.bed

            bedtools sort -i bed_files/enz_sites_sample_reads.bed > bed_files/tmp.bed

            mv bed_files/tmp.bed sample_sites_generation/enz_sites_sample_reads.bed

            bedtools cluster -d -1 -s -i sample_sites_generation/enz_sites_sample_reads.bed > {output[0]}

            Rscript split_clusters.R --out1 {output[1]} {output[0]}

        """

rule find_unique_pos:
    input:
        expand("sample_sites_generation/{sample_site_definition}", sample_site_definition=config["sample_site_definition"])[1]
    output:
        expand("multiple_bed/{find_unique_pos}", find_unique_pos=config["find_unique_pos"])[0:2],
        expand("bed_files/{find_unique_pos}", find_unique_pos=config["find_unique_pos"])[2]
    shell:
        """

            cd multiple_bed
            rm -r clusters_without_pstI/
            rm -r good_clusters/

            awk -F'\t' '{{print > $7".bed"}}' ./../{input[0]}

            Rscript ./../find_unique_positions_from_splited_clusters.R

            cd good_clusters/

            ls |  while read filename;  do cat $filename;  done >./../good_clusters.bed

            cd ..

            cd clusters_without_pstI

            ls |  while read filename;  do cat $filename;  done >./../clusters_without_pstI.bed

            cd ./../../

            cat {output[0]} {output[1]} > {output[2]}

            bedtools sort -i {output[2]} > temp.bed
            mv temp.bed {output[2]}

            bedtools cluster -d -1 -s -i bed_files/msdartseq_positions.bed > bed_files/tmp.bed
            mv bed_files/tmp.bed bed_files/msdartseq_positions.bed

            bedtools sort -i {output[2]} > temp.bed
            mv temp.bed {output[2]}

            awk 'BEGIN{{OFS="\t"}}; {{name="MS-DArT_site_"NR; print $1,$2,$3,name,$5,$6,$8}}' {output[2]} > temp.bed
            mv temp.bed {output[2]}

        """

rule trimmomatic_to_counts:
    input:
        expand("{barcodes_files}", barcodes_files=config["barcodes_files"]),
        expand("{adapters_file}", adapters_file=config["adapters_file"]),
        expand("{fastq}.FASTQ", fastq=config["fastq"])
    output:
        expand("{fastq}.without_barc_and_adapt_plus_qc.fq", fastq=config["fastq"])
    shell:
        """

            # defines the input files

            barcode_file1={input[0]}
            barcode_file2={input[1]}
            adapters_file={input[2]}

            # execute using the two sets of barcodes

            while read FILE barcode_sequence barcode_size
            do
            trimmomatic SE -threads 12 $FILE.FASTQ $FILE.without_barc_and_adapt_plus_qc.fq ILLUMINACLIP:$adapters_file:2:30:6 SLIDINGWINDOW:5:25 HEADCROP:$barcode_size MINLEN:20
            done < "$barcode_file1"

            while read FILE barcode_sequence barcode_size
            do
            trimmomatic SE -threads 12 $FILE.FASTQ $FILE.without_barc_and_adapt_plus_qc.fq ILLUMINACLIP:$adapters_file:2:30:6 SLIDINGWINDOW:5:25 HEADCROP:$barcode_size MINLEN:20
            done < "$barcode_file2"
        """

rule fastqc_to_counts:
    input:
        expand("{fastq}.without_barc_and_adapt_plus_qc.fq", fastq=config["fastq"])
    output:
        expand("fastqc_after_trimming/{fastq}.without_barc_and_adapt_plus_qc_fastqc.html", fastq=config["fastq"]),
        expand("fastqc_after_trimming/{fastq}.without_barc_and_adapt_plus_qc_fastqc.zip", fastq=config["fastq"])
    shell:
        "fastqc *_qc.fq;"
        "mv *fastqc.zip fastqc_after_trimming;"
        "mv *fastqc.html fastqc_after_trimming"

rule mapping_to_counts:
    input:
        expand("{reference_genome}", reference_genome=config["reference_genome"]),
        expand("{fastq}.without_barc_and_adapt_plus_qc.fq", fastq=config["fastq"])
    output:
        expand("mapping/{fastq}_{map}.fq", fastq=config["fastq"], map=config["mapping"]["fastq"]),
        expand("mapping/{fastq}{sam_samples}.sam", fastq=config["fastq"], sam_samples=config["mapping"]["sam_samples"])
    shell:
        """

            base=$(basename {input[0]} ".fa")

            bowtie2-build -f --threads 24 {input[0]} $base

        # 0 mismatch
        for sample in `ls ./*_qc.fq`
        do

        base=$(basename $sample ".without_barc_and_adapt_plus_qc.fq")
        echo $base
        bowtie2 -p 24 -q -k 2 -N 0 -L 15 -i S,1,0.50 -R 3 --no-1mm-upfront --end-to-end --no-unal --score-min L,0,0 -x Egrandis_297_v2.0.softmasked --un mapping/${{base}}_no_mapped_mm0.fq -U $sample -S mapping/${{base}}_mm0.sam

        done

        #1 mismatch
        for sample in `ls mapping/*_no_mapped_mm0.fq`
        do
        base=$(basename $sample "_no_mapped_mm0.fq")

        echo $base
        bowtie2 -p 24 -q -k 2 -N 0 -L 15 -i S,1,0.50 -R 3 --no-1mm-upfront --end-to-end --no-unal --mp 2 --np 10 --rdg 10 --rfg 10 --score-min L,-2,0 -x Egrandis_297_v2.0.softmasked --un mapping/${{base}}_no_mapped_mm1.fq -U mapping/${{base}}_no_mapped_mm0.fq -S mapping/${{base}}_mm1.sam
        done

        #2 mismatches
        for sample in `ls mapping/*_no_mapped_mm1.fq`
        do
        base=$(basename $sample "_no_mapped_mm1.fq")
        echo $base
        bowtie2 -p 24 -q -k 2 -N 0 -L 15 -i S,1,0.50 -R 3 --no-1mm-upfront --end-to-end --no-unal --mp 2 --np 10 --rdg 10 --rfg 10 --score-min L,-4,0 -x Egrandis_297_v2.0.softmasked --un mapping/${{base}}_no_mapped_mm2.fq -U mapping/${{base}}_no_mapped_mm1.fq -S mapping/${{base}}_mm2.sam
        done

        #3 mismatches
        for sample in `ls mapping/*_no_mapped_mm2.fq`
        do
        base=$(basename $sample "_no_mapped_mm2.fq")

        echo $base
        bowtie2 -p 24 -q -k 2 -N 0 -L 15 -i S,1,0.50 -R 3 --no-1mm-upfront --end-to-end --no-unal --mp 2 --np 10 --rdg 10 --rfg 10 --score-min L,-6,0 -x Egrandis_297_v2.0.softmasked --un mapping/${{base}}_no_mapped_mm3.fq -U mapping/${{base}}_no_mapped_mm2.fq -S mapping/${{base}}_mm3.sam
        done

        """

rule extract_sing_mapping_to_counts:
    input:
        expand("mapping/{fastq}{sam_samples}.sam", fastq=config["fastq"], sam_samples=config["mapping"]["sam_samples"])
    output:
        expand("mapping/{fastq}{sam_samples}.single_map.sam", fastq=config["fastq"], sam_samples=config["mapping"]["sam_samples"])
    shell:
        """
            for sample in `ls mapping/*.sam`
            do
            base=$(basename $sample \".sam\")

            samtools view -h mapping/${{base}}.sam | grep -v \"XS:i:\" > mapping/${{base}}.single_map.sam
            done

        """

rule sam_to_bam_to_counts:
    input:
        expand("mapping/{fastq}{sam_samples}.single_map.sam", fastq=config["fastq"], sam_samples=config["mapping"]["sam_samples"])
    output:
        expand("mapping/{fastq}{sam_samples}.single_map.bam", fastq=config["fastq"], sam_samples=config["mapping"]["sam_samples"])
    shell:
        """

            for sample in `ls mapping/*single_map.sam`
            do
            base=$(basename $sample \".sam\")
            samtools view -b mapping/${{base}}.sam > mapping/${{base}}.bam
            samtools sort mapping/${{base}}.bam > mapping/tmp.bam
            mv mapping/tmp.bam mapping/${{base}}.bam
            done

        """

rule merge_bam_to_counts:
    input:
        expand("mapping/{fastq}{sam_samples}.single_map.bam", fastq=config["fastq"], sam_samples=config["mapping"]["sam_samples"])
    output:
        expand("mapping/{fastq}_combined.bam", fastq=config["fastq"])
    shell:
        """

            for sample in `ls mapping/*_mm0.single_map.bam`
            do
            base=$(basename $sample "_mm0.single_map.bam")

            samtools merge -f mapping/${{base}}_combined.bam mapping/${{base}}_mm0.single_map.bam mapping/${{base}}_mm1.single_map.bam mapping/${{base}}_mm2.single_map.bam mapping/${{base}}_mm3.single_map.bam

            samtools sort mapping/${{base}}_combined.bam > tmp.bam
            mv tmp.bam mapping/${{base}}_combined.bam

            done

        """

rule featureCounts_input:
    input:
        expand("bed_files/{find_unique_pos}", find_unique_pos=config["find_unique_pos"])[2]
    output:
        expand("bed_files/{find_unique_pos}", find_unique_pos=config["featurecounts"])[0]
    shell:
        """
        awk 'BEGIN{{OFS="\t"}}; {{print $4,$1,$2,$3,$6}}' {input[0]} > {output[0]}
        """

rule featureCounts:
    input:
        expand("bed_files/{find_unique_pos}", find_unique_pos=config["featurecounts"])[0],
        expand("{names_correspondence}", names_correspondence=config["names_correspondence"]),
        expand("mapping/{fastq}_combined.bam", fastq=config["fastq"])
    output:
        expand("counts/{find_unique_pos}.tst", find_unique_pos=config["featurecounts"])[1]
    shell:
        """
        #rm mapping/*.sam
        #rm mapping/*single_map.sam

        declare -a list_of_bams
        list_of_bams=({input})
        list_of_bams="${{list_of_bams[@]:2}}"

        featureCounts -f -F SAF -R CORE -s 1 -T 12 -O -a {input[0]} -o {output[0]} $list_of_bams

        sed '1d' {output[0]} > tmp_counts
        mv tmp_counts {output[0]}

        while IFS=',' read bam_name sample_name
        do
        echo "sed -i 's|$bam_name|$sample_name|g' {output[0]}" >> rename.sh
        done < {input[1]}

        bash rename.sh

        rm rename.sh
        """

rule counts_correction:
    input:
        expand("bed_files/{find_unique_pos}", find_unique_pos=config["find_unique_pos"])[2],
        expand("counts/{find_unique_pos}.tst", find_unique_pos=config["featurecounts"])[1]
    output:
        expand("bed_files/{counts_correction}", counts_correction=config["counts_correction"])[0],
        expand("counts/{counts_correction}", counts_correction=config["counts_correction"])[1]
    shell:
        "Rscript counts_correction.R --out1 {output[0]} --out2 {output[1]} {input[0]} {input[1]};"

rule ALFA:
    input:
        "Egrandis_297_v2.0.gene_exons.gff3",
        "Egrandis_297_v2.0.softmasked_sizes.tst",
        expand("mapping/{fastq}_combined.bam", fastq=config["fastq"])
    output:
        "Egrandis_297_v2.0.gene_exons.gtf",
        "Egrandis_297_v2.0.gene_biotype.gtf",
        "plots_sample.Biotypes.svg",
        "plots_sample.Categories.svg"
    shell:
        """
        Rscript gff_to_gtf_biotype.R --out1 {output[0]} --out2 {output[1]} {input[0]}

        sed -i 's/UTR/utr/g' {output[1]}
        sed -i 's/mRNA/transcript/g' {output[1]}

        alfa -a {output[1]} -g Egrandis -p 7 --chr_len {input[1]}

        alfa -g ./Egrandis --bam mapping/916197_combined.bam Adult_leaves_rep1-mspI mapping/916199_combined.bam Adult_leaves_rep2-mspI mapping/916201_combined.bam Adult_leaves_rep3-mspI mapping/916198_combined.bam Juvenile_leaves_rep1-mspI mapping/916200_combined.bam Juvenile_leaves_rep2-mspI mapping/916202_combined.bam Juvenile_leaves_rep3-mspI mapping/916203_combined.bam Xylem_rep1-mspI mapping/916204_combined.bam Xylem_rep2-mspI mapping/916205_combined.bam Xylem_rep3-mspI mapping/915871_combined.bam Adult_leaves_rep1-pstII mapping/915800_combined.bam Adult_leaves_rep2-pstII mapping/915824_combined.bam Adult_leaves_rep3-pstII mapping/915882_combined.bam Juvenile_leaves_rep1-pstII mapping/915812_combined.bam Juvenile_leaves_rep2-pstII mapping/915836_combined.bam Juvenile_leaves_rep3-pstII mapping/915883_combined.bam Xylem_rep1-pstII mapping/915848_combined.bam Xylem_rep2-pstII mapping/915860_combined.bam Xylem_rep3-pstII --svg plots_samples -p 7
        """

rule marks_with_msp_bigger_than_0:
    input:
        expand("counts/{counts_correction}", counts_correction=config["counts_correction"])[1]
    output:
        expand("true_sites/{output_bigger_than_0}.txt", output_bigger_than_0=config["marks_with_msp_outputs"])
    shell:
        "Rscript marks_with_msp_bigger_than_0.R --out1 {output[0]} {input[0]}"

rule find_methylation_site_position_part1:
    input:
        expand("{met_site_input}", met_site_input=config["find_methylation_site_position_inputs_p1"])
    output:
    shell:
        "R CMD INSTALL {input[0]};"

rule find_methylation_site_position_part2:
    input:
        expand("{met_site_input}", met_site_input=config["find_methylation_site_position_inputs_p2"])
    output:
        expand("restriction_sites/{met_site_output}.bed", met_site_output=config["find_methylation_site_position_outputs_p2"])[0],
        expand("restriction_sites/{met_site_output}.bed", met_site_output=config["find_methylation_site_position_outputs_p2"])[1]
    shell:
        "Rscript restriction_sites_search.R --genome BSgenome.Egrandis.JGI.297 --out1 {output[0]} --out2 {output[1]} {input[0]}"

rule determines_sampled_site_position:
    input:
        expand("bed_files/{counts_correction}", counts_correction=config["counts_correction"])[0],
        expand("restriction_sites/{met_site_output}.bed", met_site_output=config["find_methylation_site_position_outputs_p2"])[0],
        expand("counts/{counts_correction}", counts_correction=config["counts_correction"])[1]
    output:
        expand("position_of_the_sampled_sites/{methylatio_site_output}.bed", methylatio_site_output=config["detemines_sampled_site_position"])[0],
        expand("position_of_the_sampled_sites/{methylatio_site_output}.bed", methylatio_site_output=config["detemines_sampled_site_position"])[1],
        expand("{methylatio_site_output}.csv", methylatio_site_output=config["detemines_sampled_site_position"])[2]
    shell:
        "Rscript marks_closest_restriction_site_search.R --out1 {output[0]} --out2 {output[1]} --out3 {output[2]} {input[0]} {input[1]} {input[2]}"

rule detemines_sampled_site_position_part_2:
    input:
        expand("position_of_the_sampled_sites/{methylatio_site_output}.bed", methylatio_site_output=config["detemines_sampled_site_position"])[1]
    output:
        expand("position_of_the_sampled_sites/{methylatio_site_output}_merged.tst", methylatio_site_output=config["detemines_sampled_site_position"])[1],
        expand("position_of_the_sampled_sites/{methylatio_site_output}_merged.bed", methylatio_site_output=config["detemines_sampled_site_position"])[1]
    shell:
        """

        bedtools merge -d -1 -s -c 6,2,3,4 -o distinct,collapse,collapse,collapse -delim '|' -i {input[0]} > {output[0]}

        awk 'BEGIN{{OFS=\"\t\"}}; {{print $1,$2,$3,$7,0,$4}}' {output[0]} > tempfile & mv tempfile {output[1]}

        """

rule fragments_analysis:
    input:
        expand("restriction_sites/{met_site_output}.bed", met_site_output=config["find_methylation_site_position_outputs_p2"])[0],
        expand("position_of_the_sampled_sites/{methylatio_site_output}.bed", methylatio_site_output=config["detemines_sampled_site_position"])[0]
    output:
        expand("in_silico_frags/{fragments_analysis_outptus}.bed", fragments_analysis_outptus=config["fragments_analysis_outptus"])[0],
        expand("in_silico_frags/{fragments_analysis_outptus}.txt", fragments_analysis_outptus=config["fragments_analysis_outptus"])[1],
        expand("images/in_silico_frags/{fragments_analysis_outptus}.svg", fragments_analysis_outptus=config["fragments_analysis_outptus"])[2]
    shell:
        "Rscript fragments_analysis.R --out1 {output[0]} --out2 {output[1]} --out3 {output[2]} {input[0]} {input[1]}"

rule distribution_of_the_selected_marks:
    input:
        expand("position_of_the_sampled_sites/{methylatio_site_output}_merged.bed", methylatio_site_output=config["detemines_sampled_site_position"])[1]
    output:
        expand("images/sites_distribution/{dist_graph}", dist_graph=config["dist_graph"])[0]
    shell:
        "Rscript marks_dist_chr_plot.R --out1 {output[0]} {input[0]}"

rule intersect_marks:
    input:
        expand("position_of_the_sampled_sites/{methylatio_site_output}_merged.bed", methylatio_site_output=config["detemines_sampled_site_position"])[1],
        expand("true_sites/{output_bigger_than_0}.txt", output_bigger_than_0=config["marks_with_msp_outputs"])
    params:
        names=expand("{intersect_params}", intersect_params=config["clones_names"]),
        tissue=expand("{intersect_params}", intersect_params=config["tissues_g4"]),
        enzyme=expand("{intersect_params}", intersect_params=config["enzymes"]),
        prefix=expand("{intersect_params}", intersect_params=config["intersect_marks_params"]["prefix"]),
        grupos_intersect=expand("true_sites/{output_bigger_than_0}.txt", output_bigger_than_0=config["marks_with_msp_outputs"])
    output:
        expand("{intersect_params}_intersect_marks.txt", intersect_params=config["intersect_marks_params"]["prefix"])
    script:
        "intersect_marks.R"

rule distribution_of_the_marks_in_intersection:
    input:
        expand("position_of_the_sampled_sites/{methylatio_site_output}_merged.bed", methylatio_site_output=config["detemines_sampled_site_position"])[1],
        expand("{intersect_params}_intersect_marks.txt", intersect_params=config["intersect_marks_params"]["prefix"])
    output:
        expand("images/sites_distribution/{dist_graph}", dist_graph=config["dist_graph"])[1]
    script:
        "intersect_marks_dist_chr_plot.R"

rule edgeR_with_DArTCounts:
    input:
        expand("{intersect_params}_intersect_marks.txt", intersect_params=config["intersect_marks_params"]["prefix"]),
        expand("{methylatio_site_output}.csv", methylatio_site_output=config["detemines_sampled_site_position"])[2],
    params:
        groups=expand("{edgeR_params}", edgeR_params=config["edgeR_with_DArTCounts_params"]["groups"]),
        prefix=expand("{edgeR_params}", edgeR_params=config["edgeR_with_DArTCounts_params"]["prefix"]),
        genotypes=expand("{edgeR_params}", edgeR_params=config["clones_names"]),
        tissues=expand("{edgeR_params}", edgeR_params=config["tissues_g4"]),
        sep_into=expand("{edgeR_params}", edgeR_params=config["edgeR_with_DArTCounts_params"]["sep_into"]),
        subset_model=expand("{edgeR_params}", edgeR_params=config["edgeR_with_DArTCounts_params"]["subset_model"]),
        no_bio_rep=expand("{edgeR_params}", edgeR_params=config["edgeR_with_DArTCounts_params"]["no_bio_rep"]),
        dispersion=expand("{edgeR_params}", edgeR_params=config["edgeR_with_DArTCounts_params"]["dispersion"]),
        min_msp=expand("{edgeR_params}", edgeR_params=config["edgeR_with_DArTCounts_params"]["min_msp"]),
        fdr=expand("{edgeR_params}", edgeR_params=config["edgeR_with_DArTCounts_params"]["fdr"]),
        log_fold_change=expand("{edgeR_params}", edgeR_params=config["edgeR_with_DArTCounts_params"]["log_fold_change"]),
        filtration_mode=expand("{edgeR_params}", edgeR_params=config["edgeR_with_DArTCounts_params"]["filter"]),
        number_of_tec_rep=expand("{edgeR_params}", edgeR_params=config["edgeR_with_DArTCounts_params"]["number_of_tec_rep"]),
        samples_with_tec_reps=expand("{edgeR_params}", edgeR_params=config["edgeR_with_DArTCounts_params"]["samples_with_tec_reps"]),
        samples_without_rep=expand("{edgeR_params}", edgeR_params=config["edgeR_with_DArTCounts_params"]["samples_without_rep"])
    output:
        expand("{prefix_edgeR}{groups}_DE_stats.txt", prefix_edgeR=config["edgeR_with_DArTCounts_params"]["prefix"], groups=config["edgeR_with_DArTCounts_params"]["groups"]),
        expand("{prefix_edgeR}{groups}_DE_marks.txt", prefix_edgeR=config["edgeR_with_DArTCounts_params"]["prefix"], groups=config["edgeR_with_DArTCounts_params"]["groups"]),
        expand("{prefix_edgeR}{groups}_dispersions.txt", prefix_edgeR=config["edgeR_with_DArTCounts_params"]["prefix"], groups=config["edgeR_with_DArTCounts_params"]["groups"]),
        expand("BRASUZ1_{tissue}_msp_bigger_than_threshold.txt", tissue=config["tissues_g4"])
    script:
        "edgeR_with_DArTCounts.R"

rule DEseq2_with_DArTCounts:
    input:
        expand("{intersect_params}_intersect_marks.txt", intersect_params=config["intersect_marks_params"]["prefix"]),
        expand("{methylatio_site_output}.csv", methylatio_site_output=config["detemines_sampled_site_position"])[2]
    params:
        groups=expand("{deseq_params}", deseq_params=config["DEseq2_with_DArTCounts_params"]["groups"]),
        prefix=expand("{deseq_params}", deseq_params=config["DEseq2_with_DArTCounts_params"]["prefix"]),
        genotypes=expand("{deseq_params}", deseq_params=config["clones_names"]),
        tissues=expand("{deseq_params}", deseq_params=config["tissues_g4"]),
        sep_into=expand("{deseq_params}", deseq_params=config["DEseq2_with_DArTCounts_params"]["sep_into"]),
        subset_model=expand("{deseq_params}", deseq_params=config["DEseq2_with_DArTCounts_params"]["subset_model"]),
        no_bio_rep=expand("{deseq_params}", deseq_params=config["DEseq2_with_DArTCounts_params"]["no_bio_rep"]),
        dispersion=expand("{deseq_params}", deseq_params=config["DEseq2_with_DArTCounts_params"]["dispersion"]),
        min_msp=expand("{deseq_params}", deseq_params=config["DEseq2_with_DArTCounts_params"]["min_msp"]),
        fdr=expand("{deseq_params}", deseq_params=config["DEseq2_with_DArTCounts_params"]["fdr"]),
        log_fold_change=expand("{deseq_params}", deseq_params=config["DEseq2_with_DArTCounts_params"]["log_fold_change"]),
        filtration_mode=expand("{deseq_params}", deseq_params=config["DEseq2_with_DArTCounts_params"]["filter"]),
        number_of_tec_rep=expand("{deseq_params}", deseq_params=config["DEseq2_with_DArTCounts_params"]["number_of_tec_rep"]),
        samples_with_tec_reps=expand("{deseq_params}", deseq_params=config["DEseq2_with_DArTCounts_params"]["samples_with_tec_reps"]),
        samples_without_rep=expand("{deseq_params}", deseq_params=config["DEseq2_with_DArTCounts_params"]["samples_without_rep"])
    output:
        expand("{prefix_deseq}{groups}_DE_stats.txt", prefix_deseq=config["DEseq2_with_DArTCounts_params"]["prefix"], groups=config["DEseq2_with_DArTCounts_params"]["groups"]),
        expand("{prefix_deseq}{groups}_DE_marks.txt", prefix_deseq=config["DEseq2_with_DArTCounts_params"]["prefix"], groups=config["DEseq2_with_DArTCounts_params"]["groups"])
    script:
        "deseq2_with_DArTCounts.R"

rule edger_vs_deseq:
    input:
        expand("{prefix_edgeR}{groups}_DE_marks.txt", prefix_edgeR=config["edgeR_with_DArTCounts_params"]["prefix"], groups=config["edgeR_with_DArTCounts_params"]["groups"]),
        expand("{prefix_deseq}{groups}_DE_marks.txt", prefix_deseq=config["DEseq2_with_DArTCounts_params"]["prefix"], groups=config["DEseq2_with_DArTCounts_params"]["groups"]),
        expand("position_of_the_sampled_sites/{methylatio_site_output}_merged.bed", methylatio_site_output=config["detemines_sampled_site_position"])[1]
    output:
        expand("edger_vs_deseq/{edger_vs_deseq2}", edger_vs_deseq2=config["edger_vs_deseq2_output"])[0],
        expand("edger_vs_deseq/{edger_vs_deseq2}", edger_vs_deseq2=config["edger_vs_deseq2_output"])[1],
        expand("images/edger_vs_deseq/{edger_vs_deseq2}", edger_vs_deseq2=config["edger_vs_deseq2_output"])[2]
    script:
        "edger_vs_deseq2.R"

rule counts_larger_than_threshold:
    input:
        expand("BRASUZ1_{tissue}_msp_bigger_than_threshold.txt", tissue=config["tissues_g4"])
    output:
        expand("{file}", file=config["counts_larger_than_threshold"])
    shell:
        """
        cat {input} > {output[0]}
        """

rule reproducibility_plot:
    input:
        expand("{file}", file=config["counts_larger_than_threshold"]),
        expand("counts/{counts_correction}", counts_correction=config["counts_correction"])[1],
        expand("edger_vs_deseq/{edger_vs_deseq2}", edger_vs_deseq2=config["edger_vs_deseq2_output"])[0]
    output:
        expand("images/reprodutibility/{reproducibility_plot_name}.svg", reproducibility_plot_name=config["reproducibility_plot_name"])[0],
        expand("{reproducibility_plot_name}.svg", reproducibility_plot_name=config["reproducibility_plot_name"])[1],
        expand("{reproducibility_plot_name}_sampled_sites.svg", reproducibility_plot_name=config["reproducibility_plot_name"])[1],
        expand("{reproducibility_plot_name}_methylated_sites.svg", reproducibility_plot_name=config["reproducibility_plot_name"])[1]
    params:
        prefix=expand("{reproducibility_plot_name}", reproducibility_plot_name=config["reproducibility_plot_name"])[1]
    shell:
        "Rscript reproducibility.R --output {output[0]} --correlation_output {params.prefix} {input[0]} {input[1]} {input[2]}"

rule make_bed_methylated_sites:
    input:
        expand("position_of_the_sampled_sites/{methylatio_site_output}_merged.bed", methylatio_site_output=config["detemines_sampled_site_position"])[1],
        expand("edger_vs_deseq/{edger_vs_deseq2}", edger_vs_deseq2=config["edger_vs_deseq2_output"])[1]
    output:
        expand("methylated_sites/{make_bed_methylated_sites}", make_bed_methylated_sites=config["make_bed_methylated_sites"])
    script:
        "make_bed_methylated_sites.R"

rule genomic_context:
    input:
        expand("methylated_sites/{make_bed_methylated_sites}", make_bed_methylated_sites=config["make_bed_methylated_sites"])[1],
        expand("edger_vs_deseq/{edger_vs_deseq2}", edger_vs_deseq2=config["edger_vs_deseq2_output"])[1],
        expand("{genome}", genome=config["genome_annotation_file"]),
        expand("{transposons_file}", transposons_file=config["transposons_file"])
    output:
        expand("genomic_context_files/{genomic_context}", genomic_context=config["genomic_context_outputs"])[0:6],
        expand("images/genomic_context/{genomic_context}", genomic_context=config["genomic_context_outputs"])[6],
        expand("{p2_out2}", p2_out2=config["p2_out2"]),
        expand("{p4_out1}", p4_out1=config["p4_out1"]),
        expand("{p4_out2}", p4_out2=config["p4_out2"]),
        expand("{p4_out3}", p4_out3=config["p4_out3"]),
        expand("{p4_out4}", p4_out4=config["p4_out4"])
    params:
        p1_out1=config["p1_out1"],
        p1_out2=config["p1_out2"],
        p2_out1=config["p2_out1"],
        p2_out2=config["p2_out2"],
        p3_out1=config["p3_out1"],
        p3_out2=config["p3_out2"],
        p3_out3=config["p3_out3"],
        p4_out1=config["p4_out1"],
        p4_out2=config["p4_out2"],
        p4_out3=config["p4_out3"],
        p4_out4=config["p4_out4"]
    shell:
        #Prepares the files
        "grep \"gene\" {input[2]} > genomic_context_files/gene_features.gff3;"
        "sortBed -i genomic_context_files/gene_features.gff3 > {output[0]};"
        "grep \"exon\" {input[2]} > genomic_context_files/exons_features.gff3;"
        "intersectBed -nonamecheck -wo -a {input[0]} -b {output[0]} > {output[1]};"
        #Part1
        "Rscript id_part1.R --out1 {params.p1_out1} --out2 {params.p1_out2} {input[0]} {output[1]};"
        "sortBed -i {params.p1_out2} > {params.p1_out2}.sorted.bed;"
        "intersectBed -nonamecheck -wo -a {params.p1_out2}.sorted.bed -b {input[3]} > {output[2]};"
        #Part2
        "Rscript id_part2.R --out1 {params.p2_out1} --out2 {params.p2_out2} {params.p1_out2} {output[2]};"
        "intersectBed -nonamecheck -wo -a {params.p1_out1} -b genomic_context_files/exons_features.gff3 > {output[3]};"
        #Part3
        "Rscript id_part3.R --out1 {params.p3_out1} --out2 {params.p3_out2} --out3 {params.p3_out3} {params.p1_out1} {output[3]};"
        "intersectBed -nonamecheck -wo -a {params.p3_out1} -b {input[3]} > {output[4]};"
        "intersectBed -nonamecheck -wo -a {params.p3_out3} -b {input[3]} > {output[5]};"
        #Part4
        "Rscript id_part4.R --out1 {params.p4_out1} --out2 {params.p4_out2} --out3 {params.p4_out3} --out4 {params.p4_out4} {params.p3_out1} {output[4]} {params.p3_out3} {output[5]};"
        #"Rscript id_pie_chart.R"
        "Rscript id_pie_chart.R --out1 {output[6]} {params.p2_out2} {params.p2_out1} {params.p4_out1} {params.p4_out2} {params.p3_out2} {params.p4_out3} {params.p4_out4} {input[1]}"

rule distance_to_genes:
    input:
        expand("{p2_out2}", p2_out2=config["p2_out2"]),
        expand("genomic_context_files/{genomic_context}", genomic_context=config["genomic_context_outputs"])[0]
    output:
        expand("distance_to_genes_and_TEs/{distance_to_genes_outputs}", distance_to_genes_outputs=config["distance_to_genes_outputs"])
    shell:
        "sortBed -i {input[0]} > marks_sorted.bed;"
        "closestBed -nonamecheck -s -D b -a marks_sorted.bed -b {input[1]} > {output[0]};"

rule distance_to_transposons:
    input:
        expand("{p2_out2}", p2_out2=config["p2_out2"]),
        expand("{transposons_file}", transposons_file=config["transposons_file"])
    output:
        expand("distance_to_genes_and_TEs/{distance_to_transposons_outputs}", distance_to_transposons_outputs=config["distance_to_transposons_outputs"])
    shell:
        "sortBed -i {input[0]} > marks_sorted.bed;"
        "closestBed -nonamecheck -s -D b -a marks_sorted.bed -b {input[1]} > {output[0]};"

rule restriction_sites_distance_to_genes_and_TEs:
    input:
        expand("restriction_sites/{met_site_output}.bed", met_site_output=config["find_methylation_site_position_outputs_p2"])[1],
        "Egrandis_297_v2.0_transposons.bed",
        expand("genomic_context_files/{genomic_context}", genomic_context=config["genomic_context_outputs"])[0]
    output:
        expand("distance_to_genes_and_TEs/{mspI_dist}.tst", mspI_dist=config["restriction_sites_distance_to_genes_and_TEs"])[0],
        expand("distance_to_genes_and_TEs/{mspI_dist}.tst", mspI_dist=config["restriction_sites_distance_to_genes_and_TEs"])[1]
    shell:
        "sortBed -i {input[0]} > restriction_sites/tmp.bed;"
        "mv restriction_sites/tmp.bed {input[0]};"
        "closestBed -nonamecheck -s -D b -a {input[0]} -b {input[1]} > {output[0]};"
        "closestBed -nonamecheck -s -D b -a {input[0]} -b {input[2]} > {output[1]}"

rule distance_graphs_g4:
    input:
        expand("distance_to_genes_and_TEs/{distance_to_genes_outputs}", distance_to_genes_outputs=config["distance_to_genes_outputs"])[0],
        expand("distance_to_genes_and_TEs/{distance_to_transposons_outputs}", distance_to_transposons_outputs=config["distance_to_transposons_outputs"]),
        expand("edger_vs_deseq/{edger_vs_deseq2}", edger_vs_deseq2=config["edger_vs_deseq2_output"])[1],
        expand("distance_to_genes_and_TEs/{mspI_dist}.tst", mspI_dist=config["restriction_sites_distance_to_genes_and_TEs"])[1],
        expand("distance_to_genes_and_TEs/{mspI_dist}.tst", mspI_dist=config["restriction_sites_distance_to_genes_and_TEs"])[0]
    output:
        expand("images/distance_to_genes_and_TEs/{plots_distance}.svg", plots_distance=config["distance_to_genes_and_transposon_outputs"])
    script:
        "distance_to_genes_or_transposon.R"

rule transposons_plots:
    input:
        expand("genomic_context_files/{genomic_context}", genomic_context=config["genomic_context_outputs"])[2],
        expand("genomic_context_files/{genomic_context}", genomic_context=config["genomic_context_outputs"])[5],
        expand("edger_vs_deseq/{edger_vs_deseq2}", edger_vs_deseq2=config["edger_vs_deseq2_output"])[1],
        expand("{transposons_file}", transposons_file=config["transposons_file"])
    params:
        use_intersections=expand("{use_intersections}", use_intersections=config["transposons_plots_params"]["use_intersections"]),
    output:
        expand("genomic_context_files/{transposons_for_sample}", transposons_for_sample=config["transposons_plots"]["transposons_for_sample"]),
        expand("images/transposons_plots/{transposons_venn}", transposons_venn=config["transposons_plots"]["transposons_venn"]),
        expand("genomic_context_files/{transposons_table}", transposons_table=config["transposons_plots"]["transposons_table"]),
        expand("images/transposons_plots/{transposon_bar_plot}", transposon_bar_plot=config["transposons_plots"]["transposon_bar_plot"])
    script:
        "transposons_plots.R"

rule PCAs:
    input:
        expand("{methylatio_site_output}.csv", methylatio_site_output=config["detemines_sampled_site_position"])[2],
        expand("{intersect_params}_intersect_marks.txt", intersect_params=config["intersect_marks_params"]["prefix"]),
        expand("edger_vs_deseq/{edger_vs_deseq2}", edger_vs_deseq2=config["edger_vs_deseq2_output"])[0]
    output:
        expand("images/PCAs/{outputs_pcas}", outputs_pcas=config["outputs_pcas"])
    script:
        "pcas.R"

rule venn_plots:
    input:
        expand("edger_vs_deseq/{edger_vs_deseq2}", edger_vs_deseq2=config["edger_vs_deseq2_output"])[1],
        expand("{p2_out2}", p2_out2=config["p2_out2"]),
        expand("{p2_out1}", p2_out1=config["p2_out1"]),
        expand("{p4_out3}", p4_out3=config["p4_out3"]),
        expand("{p4_out1}", p4_out1=config["p4_out1"]),
        expand("{p4_out4}", p4_out4=config["p4_out4"]),
    params:
        tissues=config["venn_plots_params"]["tissues"],
        clone_name=config["venn_plots_params"]["clone_name"],
        save_ids=config["venn_plots_params"]["save_ids"]
    output:
        expand("{outputs_venn}", outputs_venn=config["venn_plots"]["outputs"])[0],
        expand("{outputs_venn}", outputs_venn=config["venn_plots"]["outputs"])[1],
        expand("images/venn_plots/{outputs_venn}_marks.svg", outputs_venn=config["venn_plots"]["outputs"])[2],
        expand("images/venn_plots/{outputs_venn}_marks.svg", outputs_venn=config["venn_plots"]["outputs"])[3],
        expand("images/venn_plots/{outputs_venn}_marks.svg", outputs_venn=config["venn_plots"]["outputs"])[4],
        expand("images/venn_plots/{outputs_venn}_marks.svg",        outputs_venn=config["venn_plots"]["outputs"])[5],
        expand("venn_files/{outputs_venn}_marks.txt", outputs_venn=config["venn_plots"]["outputs"])[6:16],
        expand("venn_files/{outputs_venn}_marks_in_g.txt", outputs_venn=config["venn_plots"]["outputs"])[6:16]
    script:
        "venn_diag_dartcounts.R"

rule BioMart:
    input:
        expand("{genome}", genome=config["genome_annotation_file"])
    params:
        annotation_file_rda=expand("{annotation_file_rda}", annotation_file_rda=config["biomart_params"]["annotation_file_rda"]),
        its_in_listMarts=expand("{its_in_listMarts}", its_in_listMarts=config["biomart_params"]["its_in_listMarts"]),
        biomart_name=expand("{biomart_name}", biomart_name=config["biomart_params"]["biomart_name"]),
        biomart_dataset=expand("{biomart_dataset}", biomart_dataset=config["biomart_params"]["biomart_dataset"]),
        biomart_host=expand("{biomart_host}", biomart_host=config["biomart_params"]["biomart_host"]),
        biomart_vschema=expand("{biomart_vschema}", biomart_vschema=config["biomart_params"]["biomart_vschema"]),
        c_name=expand("{c_name}", c_name=config["biomart_params"]["c_name"])
    output:
        expand("annotation/biomart/{biomart_output}", biomart_output=config["biomart_output"])
    script:
        "BioMart.R"

rule sampled_mspI_fisher_test:
    input:
        expand("restriction_sites/{met_site_output}.bed", met_site_output=config["find_methylation_site_position_outputs_p2"])[0],
        expand("genomic_context_files/{genomic_context_outputs}", genomic_context_outputs=config["genomic_context_outputs"])[0],
        expand("position_of_the_sampled_sites/{methylatio_site_output}_merged.bed",
        methylatio_site_output=config["detemines_sampled_site_position"])[1],
        "tested_sites_to_methylation.txt",
        expand("methylated_sites/{make_bed_methylated_sites}", make_bed_methylated_sites=config["make_bed_methylated_sites"])
    output:
        expand("sampled_mspI_fisher_test/{sampled_mspI_fisher_test_output}", sampled_mspI_fisher_test_output=config["sampled_mspI_fisher_test_output"])
    shell:
        """
            Rscript sampled_mspI_fisher_test.R --out1 {output[0]} {input[0]}

            intersectBed -nonamecheck -wo -a {output[0]} -b {input[1]} > {output[1]}

            intersectBed -nonamecheck -wo -a {input[2]} -b {input[1]} > {output[2]}

            Rscript sampled_mspI_fisher_test_part2.R --out1 {output[3]} {input[2]} {input[3]} {output[1]} {output[2]} {input[0]} {input[4]}

        """

rule distance_methylated_sites_to_genes:
    input:
        expand("position_of_the_sampled_sites/{methylatio_site_output}_merged.bed", methylatio_site_output=config["detemines_sampled_site_position"])[1],
        expand("genomic_context_files/{genomic_context}", genomic_context=config["genomic_context_outputs"])[0],
    output:
        expand("annotation/{annotation_table_output}", annotation_table_output=config["make_annotation_table"]["annotation_table_output"])[0]
    shell:
        "closestBed -nonamecheck -s -D b -a {input[0]} -b {input[1]} > {output[0]}"

rule make_annotation_table:
    input:
        expand("{annotation_table_inputs}", annotation_table_inputs=config["make_annotation_table"]["annotation_table_inputs"]),
        expand("annotation/biomart/{biomart_output}", biomart_output=config["biomart_output"]),
        expand("annotation/{annotation_table_output}", annotation_table_output=config["make_annotation_table"]["annotation_table_output"])[0]
    output:
        expand("annotation/biomart_blast2go/{annotation_table_output}", annotation_table_output=config["make_annotation_table"]["annotation_table_output"])[2],
        expand("annotation/biomart_blast2go/{annotation_table_output}", annotation_table_output=config["make_annotation_table"]["annotation_table_output"])[1]
    script:
        "combine_annotations.R"

rule install_org_db_package:
    input:
        expand("annotation/biomart_blast2go/{annotation_table_output}", annotation_table_output=config["make_annotation_table"]["annotation_table_output"])[2]
    output:
        "org.Egrandis.eg.db"
    script:
        "makeOrgPackage.R"

rule annotation_search:
    input:
        expand("edger_vs_deseq/{edger_vs_deseq2}", edger_vs_deseq2=config["edger_vs_deseq2_output"])[1],
        expand("{p4_out1}", p4_out1=config["p4_out1"])[0],
        expand("{p4_out4}", p4_out4=config["p4_out4"])[0],
        expand("annotation/biomart_blast2go/{annotation_table_output}", annotation_table_output=config["make_annotation_table"]["annotation_table_output"])[1],
        expand("{necessary_files}", necessary_files=config["necessary_files"])[2]
    params:
        sufix_g4=config["annotation_search_params"]["sufix_g4"],
        list_by=config["annotation_search_params"]["list_by"],
        query_is=config["annotation_search_params"]["query_is"],
        query_is_file=config["annotation_search_params"]["query_is_file"],
        isoform_name=config["annotation_search_params"]["isoform_name"],
        groups=config["annotation_search_params"]["groups"],
        OrgDb_package_name=config["enrichment_params"]["OrgDb_package_name"],
        col_gene_name=config["enrichment_params"]["col_gene_name"],
        pAdjustMethod=config["enrichment_params"]["pAdjustMethod"],
        change_universe=config["enrichment_params"]["change_universe"],
        universe_subset=config["enrichment_params"]["universe_subset"],
        pvalueCutoff=config["enrichment_params"]["pvalueCutoff"],
        qvalue=config["enrichment_params"]["qvalue"],
        sim_cutoff=config["enrichment_params"]["sim_cutoff"],
        number_of_terms_by_plot=config["enrichment_params"]["number_of_terms_by_plot"],
        plot_by=config["enrichment_params"]["plot_by"]
    output:
        expand("annotation/sample_annot/{annotation_search_outputs}", annotation_search_outputs=config["annotation_search_outputs"]),
        expand("images/cluster_profiler/{annotation_search_outputs_enrich}.svg", annotation_search_outputs_enrich=config["annotation_search_enrichment_outputs"])
    script:
        "annotation_search.R"

rule venn_plots_genes:
    input:
        expand("annotation/sample_annot/{annotation_search_outputs}", annotation_search_outputs=config["annotation_search_outputs"])[0]
    params:
        tissues=config["venn_plots_params"]["tissues"],
        clone_name=config["venn_plots_params"]["clone_name"],
        save_ids=config["venn_plots_params"]["save_ids"]
    output:
        expand("images/gene_plots/{outputs_venn}_genes.svg", outputs_venn=config["venn_plots"]["outputs"])[2],
        expand("venn_files/{outputs_venn}_genes.tst", outputs_venn=config["venn_plots"]["outputs"])[6:16]
    script:
        "venn_diag_dartcounts_genes.R"

rule annotation_search_venn_sub:
    input:
        expand("annotation/biomart_blast2go/{annotation_table_output}", annotation_table_output=config["make_annotation_table"]["annotation_table_output"])[1],
        expand("venn_files/{outputs_venn}_genes.tst", outputs_venn=config["venn_plots"]["outputs"])[6:16],
        expand("{necessary_files}", necessary_files=config["necessary_files"])[2]
    params:
        sufix_g4=config["annotation_search_params"]["sufix_g4"],
        list_by=config["annotation_search_params"]["list_by"],
        query_is=config["annotation_search_params"]["query_is"],
        query_is_file=config["annotation_search_params"]["query_is_file"],
        isoform_name=config["annotation_search_params"]["isoform_name"],
        groups=config["annotation_search_params"]["groups"],
        OrgDb_package_name=config["enrichment_params"]["OrgDb_package_name"],
        col_gene_name=config["enrichment_params"]["col_gene_name"],
        pAdjustMethod=config["enrichment_params"]["pAdjustMethod"],
        change_universe=config["enrichment_params"]["change_universe"],
        universe_subset=config["enrichment_params"]["universe_subset"],
        pvalueCutoff=config["enrichment_params"]["pvalueCutoff"],
        qvalue=config["enrichment_params"]["qvalue"],
        sim_cutoff=config["enrichment_params"]["sim_cutoff"],
        number_of_terms_by_plot=config["enrichment_params"]["number_of_terms_by_plot"],
        plot_by=config["enrichment_params"]["plot_by"]
    output:
        expand("venn_files/{annotation_venn_outputs}", annotation_venn_outputs=config["annotation_venn_outputs"])
    script:
        "annotation_search_venn.R"

rule statistical_tests:
    input:
        "E_grandis_chr_sizes.bed",
        "position_of_the_sampled_sites/msdartseq_methylation_sites_of_sequenced_fragments_merged.bed",
        "sites_on_the_intersection/MS-DArT_intersect_marks.txt",
        "edger_vs_deseq/edger_DEseq2_consensus_methylated_sites.tst",
        "genomic_context_files/gene_features_sorted.gff3",
        "annotation/sample_annot/group4_methylated_genes_per_sample.tst",
        "genomic_context_files/group4_intersect_DE_transposons.txt",
        "Egrandis_297_v2.0_transposons.bed",
        "BRASUZ1 Adult Leaf_Mehylated_sites_distribution_graph_table.tst",
        "BRASUZ1 Juvenile Leaf_Mehylated_sites_distribution_graph_table.tst",
        "BRASUZ1 Xylem_Mehylated_sites_distribution_graph_table.tst",
        "distance_to_genes_and_TEs/distance_to_genes_features.txt",
        "distance_to_genes_and_TEs/distance_to_transposons_features.txt",
        "distance_to_genes_and_TEs/mspI_distance_to_genes_features.tst",
        "distance_to_genes_and_TEs/mspI_distance_to_transposons_features.tst"
    output:
        "genome_windows_1Mb.bed",
        "genome_windows_500kb.bed",
        "genome_windows_250kb.bed",
        "genome_windows_100kb.bed",
        "counts_of_genome_windows_1Mb.tst",
        "counts_of_genome_windows_500kb.tst",
        "counts_of_genome_windows_250kb.tst",
        "counts_of_genome_windows_100kb.tst",
        "images/association_plots/association_plot_TEs_with_genome.svg",
        "images/association_plots/association_plot_genomic_context.svg",
        "images/association_plots/association_plot_mspI_vicinity_genes.svg",
        "images/association_plots/association_plot_mspI_vicinity_tes.svg"
    shell:
       """
        bedtools makewindows -w 1000000 -b {input[0]} > {output[0]}
        bedtools makewindows -w 500000 -b {input[0]} > {output[1]}
        bedtools makewindows -w 250000 -b {input[0]} > {output[2]}
        bedtools makewindows -w 100000 -b {input[0]} > {output[3]}

        grep "+" {input[1]} > temp_frags.bed

        bedtools coverage -counts -a {output[0]} -b temp_frags.bed > {output[4]}
        bedtools coverage -counts -a {output[1]} -b temp_frags.bed > {output[5]}
        bedtools coverage -counts -a {output[2]} -b temp_frags.bed > {output[6]}
        bedtools coverage -counts -a {output[3]} -b temp_frags.bed > {output[7]}

        Rscript statistical_tests.R --out1 {output[8]} --out2 {output[9]} --out3 {output[10]} --out4 {output[11]} {output[4]} {output[5]} {output[6]} {output[7]} {input[2]} temp_frags.bed {input[3]} {input[4]} {input[5]} {input[6]} {input[7]} "{input[8]}" "{input[9]}" "{input[10]}" {input[11]} {input[12]} {input[13]} {input[14]}
        """

rule fisher_test:
    input:
        expand("edger_vs_deseq/{edger_vs_deseq2}", edger_vs_deseq2=config["edger_vs_deseq2_output"])[0],
        expand("{outputs_venn}", outputs_venn=config["venn_plots"]["outputs"])[0],
        expand("annotation/sample_annot/{annotation_search_outputs}", annotation_search_outputs=config["annotation_search_outputs"])[0],
        expand("{outputs_venn}", outputs_venn=config["venn_plots"]["outputs"])[1]
    output:
    script:
        "fisher_test_brasuz1.R"

#Check the methylations in genes with GbM accordingly with the classification generated with the Brasuz1 v.1.0 Niederhuth et al., 2016 (https://doi.org/10.1186/s13059-016-1059-0)
rule genes_with_GBM:
    input:
        expand("sampled_mspI_fisher_test/{sampled_mspI_fisher_test_output}", sampled_mspI_fisher_test_output=config["sampled_mspI_fisher_test_output"])[3],
        expand("{gbm_list}", gbm_list=config["genes_with_GBM"]["gbm_list"]),
        expand("sampled_mspI_fisher_test/{sampled_mspI_fisher_test_output}", sampled_mspI_fisher_test_output=config["sampled_mspI_fisher_test_output"])[1],
        expand("{intersect_params}_intersect_marks.txt", intersect_params=config["intersect_marks_params"]["prefix"]),
        expand("position_of_the_sampled_sites/{methylatio_site_output}_merged.bed", methylatio_site_output=config["detemines_sampled_site_position"])[1],
        expand("annotation/sample_annot/{annotation_search_outputs}", annotation_search_outputs=config["annotation_search_outputs"])[1:4]
    shell:
        "echo genes_with_GBM;"
        "Rscript genes_with_GBM.R {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} {input[6]} {input[7]}"

rule bs_binomial_test:
    input:
        expand("{bissulfite_bb}.bed", bissulfite_bb=config["bissulfite_bb"])
    output:
        expand("bs_seq_files/{bs_binomial_test}.bed", bs_binomial_test=config["bs_binomial_test"])
    script:
        "binomial_teste.R"

rule ms_dart_seq_validation:
    input:
        expand("position_of_the_sampled_sites/{methylatio_site_output}_merged.bed", methylatio_site_output=config["detemines_sampled_site_position"])[1],
        expand("bs_seq_files/{bs_binomial_test}.bed", bs_binomial_test=config["bs_binomial_test"]),
        expand("edger_vs_deseq/{edger_vs_deseq2}", edger_vs_deseq2=config["edger_vs_deseq2_output"])[1]
    params:
        cpgs_validation_methy_info_CpG=expand("{cpgs_validation_methy_info_CpG}", cpgs_validation_methy_info_CpG=config["cpgs_validation"])
    output:
        expand("{ms_dart_seq_validation}.bed", ms_dart_seq_validation=config["ms_dart_seq_validation"])[0],
        expand("msdart_validation/{ms_dart_seq_validation}.tst", ms_dart_seq_validation=config["ms_dart_seq_validation"])[1:3]
    shell:
        "Rscript recover_the_restriction_site_of_each_dart_site.R --out1 {output[0]} {input[0]};"
        "mkdir temp;"
        "cat {output[0]} {input[1]} > temp/clusters_dart_bs_cpgs.bed;"
        "bedtools sort -i temp/clusters_dart_bs_cpgs.bed > temp/clusters_dart_bs_cpgs_sorted.bed;"
        "clusterBed -d -1 -i temp/clusters_dart_bs_cpgs_sorted.bed > temp/msdart_bsseq_CpG_clusters.bed;"
        "Rscript cpgs_validation.R --cpg_file {params.cpgs_validation_methy_info_CpG} --out1 {output[1]} --out2 {output[2]} temp/msdart_bsseq_CpG_clusters.bed {input[2]};"
        "rm -r temp/"
