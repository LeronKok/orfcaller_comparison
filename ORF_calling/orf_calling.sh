### Index/reference preparation - only once
##Bowtie index
bowtie2-build -f \
    --threads 4 \
    --seed 24 \
    "<contaminant_fasta>" \
    "<output_index_file>"

##STAR index
STAR --runMode genomeGenerate \
    --runThreadN 12 \
    --genomeDir "<index_save_directory>" \
    --genomeFastaFiles "<genome_fasta>" \
    --genomeSAindexNbases 10
    --sjdbOverhang 29 \
    --sjdbGTFfile "<gtf>"

##Rannot annotation - in R 4.1.2
```prepare_annotation_files(annotation_directory = <save_directory>, 
                         twobit_file = <genome_twobit_format>, 
                         gtf_file = <gtf>, 
                         annotation_name = <annot_name>, 
                         forge_BSgenome = TRUE)

```

##Gedi genome reference (for PRICE)
gedi -e IndexGenome -s "<genome_fasta>" -a "<gtf>" \
    -n "<save_name>" -f "<save_directory>" \
	-o "<save_name>/<save_directory>.oml" \
    -nobowtie -nostar -nokallisto

##Ribotricer index 
python3 ribotricer prepare-orfs --gtf "<gtf>" --fasta "<genome_fasta>" \
		 --min_orf_length 9  \
        --longest --start_codons ATA,ATC,ATG,ATT,AAG,ACG,AGG,ATG,CTG,GTG,TTG \
        --prefix "<ribotricer_index_prefix>"

### Data preprocessing
##Trimgalore
trim_galore "<fastq_input>" \
  --cores 2 \
  --gzip \
  --length 25 \
  --trim-n \
  --fastqc \
  --fastqc_args "--outdir <output_directory>" \
  --output_dir "<output_directory>"

##Contaminant filtering
bowtie2 --seedlen=25 \
  --threads 12 \
  --time \
  --un-gz "<output_file>" \
  -x ${bowtie2_index} \
  -U "<trimgalore_result_file>" \
  -S "<contaminants_found.sam>"

##Mapping
STAR --genomeDir "<star_29nt_index>" \
  --sjdbGTFfile "<gtf>" \
  --runThreadN 24 \
  --runDirPerm All_RWX \
  --twopassMode Basic \
  --readFilesIn \
  "<bowtie2_filtered.fastq.gz>" \
  --readFilesCommand zcat \
  --outFilterMismatchNmax 2 \
  --outFilterMultimapNmax 20 \
  --outSAMattributes All \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts \
  --outFileNamePrefix "<star_result_prefix>" \
  --limitOutSJcollapsed 10000000 \
  --limitIObufferSize=300000000 \
  --outFilterType BySJout \
  --alignSJoverhangMin 1000 \
  --outTmpKeep None

samtools index -@ 24 "<star_output>"

### ORF calling
##ORFquant - in R 4.1.2
```RiboseQC_analysis(annotation_file = <rannot>,
                        bam_files = <star_bamfile>,
                        read_subset = F,
                        readlength_choice_method = "max_coverage",
                        dest_names = <output_dir>,
                        rescue_all_rls = FALSE, fast_mode = F,
                        create_report = T, sample_names = NA,
                        report_file = <output_dir>, extended_report = F,
                        pdf_plots = T)
```

```run_ORFquant(for_ORFquant_file = <for_ORFquant_file>,
                     annotation_file = <rannot>,
                     n_cores = 24,
                     prefix = <ORFquant_out>,
                     gene_name = NA,
                     gene_id = NA,
                     genomic_region = NA,
                     write_temp_files = T,
                     write_GTF_file = T,
                     write_protein_fasta = T,
                     interactive = T,
                     stn.orf_find.all_starts = T,
                     stn.orf_find.nostarts = F,
                     stn.orf_find.start_sel_cutoff = NA,
                     stn.orf_find.start_sel_cutoff_ave = 0.5,
                     stn.orf_find.cutoff_fr_ave = 0.5,
                     stn.orf_quant.cutoff_cums = NA,
                     stn.orf_quant.cutoff_pct = 2,
                     stn.orf_quant.cutoff_P_sites = NA,
                     unique_reads_only = F,
                     canonical_start_only = T,
                     stn.orf_quant.scaling = "total_Psites"
                     )
```

##PRICE
samtools view -b -q 5 "<star_bamfile>" | samtools sort > "<filtered_bamfile>"
samtools index "<filtered_bamfile>"

gedi -e Price -reads "<filtered_bamfile>" \
    -genomic "<gedi_genome_reference>" \
    -prefix "<output_prefix" -plot

gedi Nashorn -e 'load("'<price_orfs.cit>'").ei().map(function(o) new BedEntry(o.data.getStartStop(o,true).toMutable().setData(new NameAnnotation(o.data.getGeneId()+"__"+o.data.getTranscript()+"__"+o.data.getType()+"__"+o.data.getOrfid()+"__"+o.data.getStartCodon())))).print()' > "<price_orfs.cit.bed>"

#Extracting activity values/p-sites
gedi -e ViewCIT -m bedgraph -filter 'ref.isPlus() && d.sum()>0' -score 'd.sum()' <price_codons.cit> > <p_sites.+.bedgraph>
gedi -e ViewCIT -m bedgraph -filter 'ref.isMinus() && d.sum()>0' -score 'd.sum()' <price_codons.cit> > <p_sites.-.bedgraph>
python3 prc_process_bedgraph.sh <p_sites.+.bedgraph> <output_bedgraph> +
python3 prc_process_bedgraph.sh <p_sites.-.bedgraph> <output_bedgraph> -

##Ribo-TISH
python2 ribotish predict -b "<star_bamfile>" \
        -g "<gtf>" -f "<genome_fasta>" --longest -o "<output_file>"

##Ribotricer
python3 ribotricer detect-orfs --bam "<star_bamfile>" \
        --ribotricer_index "<ribotricer_index>" \
		--prefix "<ribotricer_results_prefix>" \
        --phase_score_cutoff 0.440