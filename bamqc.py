#!/usr/bin/env python3
__author__ = 'zoukeke'
import numpy as np
import pysam
import collections
import pandas as pd
import re
import os
import argparse

class BamQc:
    """需要排序后的bam,计算duplication会用到"""
    def __init__(self,bam_file=None,bed_file=None,ref_version=None,cond=None,outdir=os.getcwd(),is_draw_all='yes',is_draw_target='yes'):
        """
        :param bam_file:
        :param bed_file:
        :param ref_version:
        :param cond: 排序后bam文件,统计靶向区域的条件，默认为大于深度的平均值
        """
        self.bam_file = bam_file
        self.bed_file = bed_file
        self.is_draw_target = is_draw_target
        self.is_draw_all = is_draw_all
        if not ref_version:
            self.ref_len = 2897310462
        elif ref_version.lower() == "hg19" or ref_version.lower() == "grch37":
            self.ref_len = 2897310462
        elif ref_version.lower()  == "hg38" or ref_version.lower() == "grch38" :
            self.ref_len = 3049315783
        if not cond:
            self.cond = "np.mean(self.pos_depth)" # "100" "0" "np.mean(self.pos_depth)*0.5" etc. target区域位点的筛选条件
        else:
            self.cond = cond
        self.sample = os.path.basename(bam_file).split(".")[0]
        self.outdir = os.path.join(outdir,self.sample,'qc')
        os.makedirs(self.outdir,exist_ok=True)

    @staticmethod
    def __is_bedfile(bed_file):
        "通过开始和结尾是否相等判断是否为正规bed文件"
        with open(bed_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith("#"):continue
                line = line.split(None)
                chrom,start,end = line[0:3]
                if start == end :
                    return False
        return True

    def get_target(self):
        "get target region from bed file"
        target_len = 0
        if not self.bed_file:
            target = [['chr' + str(i),None,None] for i in range(1, 23)] + [['chrX',None,None], ['chrY',None,None], ['chrM',None,None]]
        else:
            target = []
            with open(self.bed_file) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("#"):continue
                    l = line.split(None)
                    chrom,start,end = l[0:3]
                    if not self.__is_bedfile(self.bed_file):
                        start = int(start) - 1
                    end = int(end)
                    target.append([chrom,int(start),int(end)])
                    target_len += int(end) - int(start)
        self.target = target
        if not self.bed_file:
            self.target_len = self.ref_len
        else:
            self.target_len = target_len


    def get_coverage(self):
        #if no target region start and stop should be None
        samfile = pysam.AlignmentFile(self.bam_file)
        self.get_target()
        mapped_alignments = samfile.mapped #-F 0x4 is_unmapped 之和
        unmapped_alignments = samfile.unmapped # -f 0x4
        self.all_alignments = mapped_alignments + unmapped_alignments #0x4

        pos_depth = np.array([])
        region_median_depth = np.array([])
        region_read = np.array([])
        x_label = []

        for chrom,start,stop in self.target: #如果是non target, 则start 和 stop都为 None
            x_label.append(chrom + ":"+str(start) + "-" + str(stop))
            x = samfile.count_coverage(chrom,start=start,stop=stop, quality_threshold=0)

            xdata = pd.DataFrame(data={'A': x[0], 'C': x[1], 'G': x[2], 'T': x[3]})
            xdata = xdata.assign(depth = lambda x: x.apply(lambda y : sum(y),axis=1))
            depth = np.array(xdata.depth)

            #cover_len += len(depth[depth>0]) #target区深度大于0的位点的个数,除以target 区域长度可以计算覆盖度
            #cover_num += np.sum(depth[depth>0]) #target区深度大于0的位点的碱基数,除以总碱基数，可以表示捕获的特异性
            pos_depth = np.append(pos_depth,depth)
            region_read = np.append(region_read,samfile.count(chrom, start=start, stop=stop))
            region_median_depth = np.append(region_median_depth,np.median(depth))


        self.cover_rate = len(pos_depth[pos_depth>0])/self.target_len #target区域大于零的覆盖度
        #self.cover_cond_rate = len(pos_depth[pos_depth > eval(self.cond)])/self.target_len #target区域大于cond的覆盖度
        self.target_mean_depth = np.mean(pos_depth)
        self.pos_depth = pos_depth
        self.region_read = region_read
        self.region_median_depth = region_median_depth
        self.x_label = x_label
        samfile.close()

    @staticmethod
    def __is_dup(last_alignment,this_alignment):
        "identify whether duplicated reads"
        if last_alignment.reference_name == this_alignment.reference_name:
            if last_alignment.reference_start == this_alignment.reference_start:
                if last_alignment.cigarstring == this_alignment.cigarstring:
                    return True
        else:
            return False

    @staticmethod
    def __get_pos(sequence=None,pos_dict=None):
        if pos_dict is None:
            pos_dict = {}
        elif not isinstance(pos_dict,dict):
            raise ValueError('must dict')
        pos = 0
        for base in sequence:
            pos += 1
            pos_dict[pos] = pos_dict.get(pos,{})
            pos_dict[pos][base] = pos_dict[pos].get(base,0)
            pos_dict[pos][base]+= 1
        return pos_dict

    @staticmethod
    def __get_q(pos_q_dict=None,v=None):
        all = 0
        for pos,q_dict in pos_q_dict.items():
            all += sum([ct for q, ct in q_dict.items() if q >= v])
        return all

    @staticmethod
    def __get_atgc(pos_atcg_dict=None):
        all = 0
        for pos,atcg_dict in pos_atcg_dict.items():
            all += sum([ct for base,ct in atcg_dict.items() if base == "G" or base == "C"])
        return all

    @staticmethod
    def __convert(odi):
        x = []
        for q, count in odi.items():
            x = x + [q] * count
        return x

    def draw_boxplot(self,od=None,pic_ty='png'):
        import matplotlib.pyplot as plt
        plt.switch_backend('agg')
        data = []
        for pos, odi in od.items():
            data.append(self.__convert(odi))

        widths = 0.1
        gap = 0.1
        positions = [num * (gap + (widths / 2)) for num in range(1, len(data) + 1)]
        fs = positions[-1] + 1
        fontsize = 1
        f, ax = plt.subplots(facecolor="w", edgecolor="w",figsize=(fs*5,fs))
        ax.set_title('quality distribution')
        # whis指定上下须与上下四分位的距离，默认为1.5倍的四分位差；
        # medianprops = dict(linestyle='-', linewidth=0.1, color='red')
        #boxprops = dict(linewidth=widths, color='black')
        flierprops = dict(marker='+', markerfacecolor='green', markersize=0.5, linestyle='none')
        ax.boxplot(data, widths=widths, showmeans=True, flierprops=flierprops, positions=positions,labels=od.keys())
        ax.tick_params(axis='both', labelsize=fontsize)
        pic = os.path.join(self.outdir,self.sample+"_quality."+pic_ty)
        f.savefig(pic, bbox_inches='tight', format=pic_ty)
        plt.cla()
        plt.close("all")


    def bam_qc(self):
        samfile = pysam.AlignmentFile(self.bam_file)

        mapped_alignments = samfile.mapped #-F 0x4 is_unmapped 之和
        unmapped_alignments = samfile.unmapped # -f 0x4
        all_alignments = mapped_alignments + unmapped_alignments #0x4

        mapping_quality = 0
        duplicate_alignments = 0   #optical or PCR duplicate, flag 1024
        PEMap = 0
        UniqueMap = 0
        insert_arr = []  #template_length list  to draw insert size plot ++
        mismatch = 0   #the sum of XM tag
        pos_atcg_dict = collections.OrderedDict()
        pos_atcg_align_dict = collections.OrderedDict()
        pos_atcg_uniq_dict = collections.OrderedDict()
        pos_quality_dict = collections.OrderedDict()

        I_cigar = 1; D_cigar = 2
        insert = 0
        deletion = 0
        gap = 0
        read_len = [1000,0]

        last_alignment = None
        for alignment in samfile:

            sequence = alignment.seq
            ln = len(sequence)
            if ln > read_len[1]: read_len[1] = ln
            if ln < read_len[0]: read_len[0] = ln
            pos_atcg_dict = self.__get_pos(sequence=sequence,pos_dict=pos_atcg_dict)
            pos_quality_dict = self.__get_pos(sequence=alignment.get_forward_qualities(),pos_dict=pos_quality_dict)
            if alignment.has_tag('XM'):
                mismatch += alignment.get_tag('XM')

            if alignment.is_unmapped:
                last_alignment = alignment
                continue

            if not last_alignment: #segment mapping 上的duplication
                last_alignment = alignment
            else:
                if self.__is_dup(last_alignment, alignment):
                    duplicate_alignments += 1
            #gc_align_dict = self.__get_gc(sequence=sequence,gc_dict=gc_align_dict)
            pos_atcg_align_dict = self.__get_pos(sequence=sequence,pos_dict=pos_atcg_align_dict)

            insert += sum([i[1] for i in alignment.cigartuples if i[0] == I_cigar])
            deletion += sum([i[1] for i in alignment.cigartuples if i[0] == D_cigar])
            if alignment.has_tag('XG'):
                gap += alignment.get_tag('XG')

            if not alignment.mate_is_unmapped:
                PEMap += 1
                if not alignment.is_secondary:
                    UniqueMap += 1
                    #gc_uniq_dict = self.__get_gc(sequence=sequence,gc_dict=gc_uniq_dict)
                    pos_atcg_uniq_dict = self.__get_pos(sequence=sequence, pos_dict=pos_atcg_uniq_dict)

            mapping_quality += alignment.mapping_quality
            insert_arr.append(abs(alignment.template_length)) #the observed query template length
            last_alignment = alignment    # final
        samfile.close()

        # result values
        mapping_rate = mapped_alignments / all_alignments  ##
        mapping_quality = mapping_quality / all_alignments ##
        duplication_rate = duplicate_alignments / all_alignments ##

        insert_arr = np.array(insert_arr)
        insert_arr = insert_arr[insert_arr>0]
        insert_size = int(np.median(insert_arr)) ##template_length > 0的中位数

        self.mapping_rate = mapping_rate
        self.mapping_quality = mapping_quality
        self.duplication_rate = duplication_rate
        self.insert_size = insert_size
        self.reads = all_alignments
        self.PEMap_rate =  PEMap / all_alignments
        self.UniqueMap_rate = UniqueMap / all_alignments
        self.bases = self.__get_q(pos_q_dict=pos_quality_dict, v=0)
        self.Q20_rate = self.__get_q(pos_q_dict=pos_quality_dict, v=20)/self.bases
        self.Q30_rate = self.__get_q(pos_q_dict=pos_quality_dict, v=30)/self.bases
        self.gc_rate = self.__get_atgc(pos_atcg_dict=pos_atcg_dict)/self.bases
        self.gc_align_rate = self.__get_atgc(pos_atcg_dict=pos_atcg_align_dict)/self.bases
        self.gc_uniq_rate = self.__get_atgc(pos_atcg_dict=pos_atcg_uniq_dict)/self.bases
        self.mismatch_rate = mismatch/self.bases #XM
        self.deletion_rate = deletion/self.bases #cigar D
        self.insert_rate = insert/self.bases #cigar I
        self.gap_rate = gap/self.bases #XG tag
        self.read_len  = read_len
        #draw picture
        self.pos_quality_dict = pos_quality_dict
        self.pos_atcg_dict = pos_atcg_dict
        self.pos_atcg_uniq_dict = pos_atcg_uniq_dict
        self.pos_atcg_align_dict = pos_atcg_align_dict
        self.insert_arr = insert_arr

    def draw_hist(self,od=None,pic_ty="png"):
        import seaborn as sns
        import matplotlib.pyplot as plt
        from scipy.stats import norm
        plt.switch_backend('agg')
        od = od[(od>0) & (od<=600)]
        f, all_pic = plt.subplots()

        pp=sns.distplot(od, fit=norm, rug=True, color="#8B4513")
        pp.set(xlim=(0, 600))
        all_pic.set_title('the distribution of insert size', fontsize=15)
        pic = os.path.join(self.outdir, self.sample + "_insert." + pic_ty)
        f.savefig(pic, bbox_inches='tight', format=pic_ty)
        plt.cla()
        plt.close("all")

    @staticmethod
    def __split(key):
        chrom,start,end = re.split(":|-",key)
        chrom = chrom.replace("chr",'')
        if chrom == 'X':
            chrom = 23
        elif chrom == 'Y':
            chrom = 24
        elif chrom == 'M':
            chrom = 0
        else:
            chrom = int(chrom)
        return chrom,int(start),int(end)

    def draw_target(self,L=None,pic_ty='png'): #L 一个target区域一个值
        import matplotlib.pyplot as plt
        import seaborn as sns
        plt.switch_backend('agg')
        step = 20
        split_pos = list(map(self.__split,self.x_label))
        chroms = [chrom for chrom,start,end in split_pos]
        starts = [start for chrom,start,end in split_pos]
        ends = [end for chrom, start, end in split_pos]
        df = (pd.DataFrame(data={"num": L,"chrom":chroms,"start":starts,'end':ends},index=self.x_label)
                .sort_values(by=['chrom','start'])
              )
        add = step - df.shape[0] % step
        num = np.array(list(df.num) + [0] * add).astype('int32')
        num = num.reshape(-1,step)
        name = np.array(list(df.index) + ['0'] * add)
        name = name.reshape(-1,step)

        f, all_pic = plt.subplots(figsize=(22, 0.5 * len(num)), nrows=1)
        #sns.heatmap(num, linewidths=.5, cmap="YlGnBu", ax=all_pic,cbar=False,annot=name,fmt="s",annot_kws={'size': 4})
        sns.heatmap(num, linewidths=.5, cmap="YlGnBu", ax=all_pic, cbar=False, annot=True, fmt=".2g")#,annot_kws={'size': 5})
        all_pic.set_title('target region coverage', fontsize=10)
        # for item in all_pic.get_xticklabels():
        #     item.set_rotation(0)
        pic = os.path.join(self.outdir, self.sample + "_target." + pic_ty)
        f.savefig(pic, bbox_inches='tight', format=pic_ty)
        plt.cla()
        plt.close("all")

    def draw_gc(self,od=None,pic_ty='png'):
        import collections
        import numpy as np
        import seaborn as sns
        import matplotlib.pyplot as plt
        import pandas as pd
        plt.switch_backend('agg')
        sns.set(style="darkgrid")
        def convert_dataframe(atgc=None):
            pos = list(atgc.keys())
            A = [base_dic['A'] / sum(base_dic.values()) for pos, base_dic in atgc.items()]
            T = [base_dic['T'] / sum(base_dic.values()) for pos, base_dic in atgc.items()]
            C = [base_dic['C'] / sum(base_dic.values()) for pos, base_dic in atgc.items()]
            G = [base_dic['G'] / sum(base_dic.values()) for pos, base_dic in atgc.items()]
            dfA = pd.DataFrame(data={'pos': pos, "count": A}).assign(base=lambda x: np.array('A'))
            dfT = pd.DataFrame(data={'pos': pos, "count": T}).assign(base=lambda x: np.array('T'))
            dfG = pd.DataFrame(data={'pos': pos, "count": G}).assign(base=lambda x: np.array('G'))
            dfC = pd.DataFrame(data={'pos': pos, "count": C}).assign(base=lambda x: np.array('C'))
            df = dfA.append(dfT, ignore_index=True)
            df = df.append(dfC, ignore_index=True)
            df = df.append(dfG, ignore_index=True)
            return df

        f, all_pic = plt.subplots()
        df = convert_dataframe(od)
        sns.lineplot(x='pos', y='count', hue="base", style="base", data=df, ax=all_pic)
        all_pic.set_title('position ATCG percent', fontsize=15)
        pic = os.path.join(self.outdir, self.sample + "_gc." + pic_ty)
        f.savefig(pic, bbox_inches='tight', format=pic_ty)
        plt.cla()
        plt.close("all")

    def get_target_stat(self):
        ##target region
        self.get_target()
        self.get_coverage()

        # self.target_len #靶向区域的长度，若无bed文件提供则为整个参考基因组区域
        # self.cover_rate #target区域大于零的覆盖度
        # self.target_mean_depth #target区域的平均覆盖深度
        # self.pos_depth #array , target靶向区域的每个pos的深度
        # self.region_read #array, 覆盖在target区域的reads(一个区域的总reads), ++ 画一个可视化的图
        # self.region_median_depth #array, 每个region上位点bases数的中位数，++ 画一个可视化的图
        self.cover_cond_rate = len(self.pos_depth[self.pos_depth > eval(self.cond)])/self.target_len #target区域大于条件的覆盖度，默认为平均值
        self.cover_cond_rate_20 = len(self.pos_depth[self.pos_depth > eval(self.cond) * 0.2]) / self.target_len  # target区域大于条件的20%的覆盖度，默认为平均值
        self.cover_cond_rate_50 = len(self.pos_depth[self.pos_depth > eval(self.cond) * 0.5]) / self.target_len  # target区域大于条件的20%的覆盖度，默认为平均值
        self.cover_30_rate = len(self.pos_depth[self.pos_depth > 30]) /self.target_len
        self.cover_100_rate = len(self.pos_depth[self.pos_depth > 100]) / self.target_len
        self.capture_bases = np.sum(self.pos_depth)
        self.capture_reads = np.sum(self.region_read)
        self.capture_rate_bases = self.capture_bases/self.all_alignments
        self.capture_rate_reads = self.capture_reads /self.all_alignments
        if self.is_draw_target == "yes":
            #self.draw_target(L=self.region_read,pic_ty='png') #picture
            self.draw_target(L=self.region_median_depth, pic_ty="png") #picture

    def get_all_stat(self):
        self.bam_qc()
        if self.is_draw_all == "yes":
            self.draw_gc(od=self.pos_atcg_dict,pic_ty='png')
            self.draw_boxplot(od=self.pos_quality_dict,pic_ty='png')
            self.draw_hist(od=self.insert_arr,pic_ty='png')

    def get_stat(self):
        self.get_all_stat()
        with open(os.path.join(self.outdir,"stat.txt"),'w') as f:
            print("reads",self.reads,sep="\t",file=f)
            print("bases",self.bases,sep="\t",file=f)
            print("len_reads",self.read_len,sep="\t",file=f)
            print("insert_mid",self.insert_size,sep="\t",file=f)
            print("Q20_rate",self.Q20_rate,sep="\t",file=f)
            print("Q30_rate",self.Q30_rate,sep="\t",file=f)
            print("gc_rate",self.gc_rate,sep="\t",file=f)
            print("gc_align_rate", self.gc_align_rate, sep="\t", file=f)
            print("gc_uniq_rate", self.gc_uniq_rate, sep="\t", file=f)
            print("mismatch_rate", self.mismatch_rate, sep="\t", file=f)
            print("deletion_rate", self.deletion_rate, sep="\t", file=f)
            print("insert_rate", self.insert_rate, sep="\t", file=f)
            print("gap_rate", self.gap_rate, sep="\t", file=f)
            print("map_quali", self.mapping_quality, sep="\t", file=f)
            print("mapping_rate", self.mapping_rate, sep="\t", file=f)
            print("UniqueMap_rate", self.UniqueMap_rate, sep="\t", file=f)
            print("duplication_rate", self.duplication_rate, sep="\t", file=f)
            print("PEMap_rate", self.PEMap_rate, sep="\t", file=f)
            print("insert_rate", self.insert_rate, sep="\t", file=f)

        if self.bed_file:
            self.get_target_stat()
            with open(os.path.join(self.outdir,"stat.txt"),'a') as f:
                print("target_len", self.target_len, sep="\t", file=f)
                print("cover_rate", self.cover_rate, sep="\t", file=f)
                print("target_mean_depth", self.target_mean_depth, sep="\t", file=f)
                print("cover_cond_rate_20", self.cover_cond_rate_20, sep="\t", file=f)
                print("cover_cond_rate_50", self.cover_cond_rate_50, sep="\t", file=f)
                print("cover_30X_rate",self.cover_30_rate,sep="\t",file=f)
                print("cover_100X_rate",self.cover_100_rate,sep="\t",file=f)
                print("capture_rate_bases", self.capture_rate_bases, np.sum(self.pos_depth), sep="\t", file=f)
                print("capture_rate_reads", self.capture_rate_reads, np.sum(self.region_read), sep="\t", file=f)

def get_par():
    #self,bam_file=None,bed_file=None,ref_version=None,cond=None,outdir=os.getcwd()
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bamfile", help="sorted bam file", default=None,required=True)
    parser.add_argument("-r", "--bedfile", help="targeted bed file", default=None)
    parser.add_argument("--rv",dest="ref_version",default="hg19",help="your reference type, default hg19")
    parser.add_argument("--cond",default="np.mean(self.pos_depth)",help="condition of cover_cond_rate, default is 'np.mean(self.pos_depth)'")
    parser.add_argument("--outdir",default=os.getcwd(),help="output directory, default is current directory")
    parser.add_argument("--is_draw_target",default="yes",help="whether draw target region pic,defult yes")
    parser.add_argument("--is_draw_all",default="yes",help="whether draw all bam pic,default yes")
    parser.add_argument("--tr",default="no",help="whether compute target report,default no")
    parser.add_argument("--ar",default="no",help="whether compute all report,default no")
    args = parser.parse_args()
    return args

def main():
    par = get_par()
    oj = BamQc(bam_file=par.bamfile,bed_file=par.bedfile,cond=par.cond,ref_version=par.ref_version,outdir=par.outdir,is_draw_target=par.is_draw_target,is_draw_all=par.is_draw_all)
    #oj.get_stat()
    if par.tr == "yes":
        m = 'w'
        mm = 'a'
    else:
        mm = 'w'
    if par.tr == "yes":
        oj.get_target_stat()
        with open(os.path.join(oj.outdir,"targetstat.txt"),m) as f:
            print("target_len", oj.target_len, sep="\t", file=f) 
            print("cover_rate", oj.cover_rate, sep="\t", file=f)
            print("target_mean_depth", oj.target_mean_depth, sep="\t", file=f)                           
            print("cover_cond_rate_20", oj.cover_cond_rate_20, sep="\t", file=f)                         
            print("cover_cond_rate_50", oj.cover_cond_rate_50, sep="\t", file=f)
            print("cover_30X_rate",oj.cover_30_rate,sep="\t",file=f)
            print("cover_100X_rate",oj.cover_100_rate,sep="\t",file=f)
            print("capture_rate_bases", oj.capture_rate_bases, sep="\t", file=f)
            print("capture_rate_reads", oj.capture_rate_reads, sep="\t", file=f)
            print("capture_bases", int(oj.capture_bases), sep="\t", file=f)
            print("capture_reads", int(oj.capture_reads), sep="\t", file=f)
    if par.ar == "yes":
        oj.get_all_stat()
        with open(os.path.join(oj.outdir,"stat.txt"),mm) as f:
            print("reads",oj.reads,sep="\t",file=f)
            print("bases",oj.bases,sep="\t",file=f)
            print("len_reads",oj.read_len,sep="\t",file=f)
            print("insert_mid",oj.insert_size,sep="\t",file=f)
            print("Q20_rate",oj.Q20_rate,sep="\t",file=f)
            print("Q30_rate",oj.Q30_rate,sep="\t",file=f)
            print("gc_rate",oj.gc_rate,sep="\t",file=f)
            print("gc_align_rate", oj.gc_align_rate, sep="\t", file=f)
            print("gc_uniq_rate", oj.gc_uniq_rate, sep="\t", file=f)
            print("mismatch_rate", oj.mismatch_rate, sep="\t", file=f)
            print("deletion_rate", oj.deletion_rate, sep="\t", file=f)
            print("insert_rate", oj.insert_rate, sep="\t", file=f)
            print("gap_rate", oj.gap_rate, sep="\t", file=f)
            print("map_quali", oj.mapping_quality, sep="\t", file=f)
            print("mapping_rate", oj.mapping_rate, sep="\t", file=f)
            print("UniqueMap_rate", oj.UniqueMap_rate, sep="\t", file=f)
            print("duplication_rate", oj.duplication_rate, sep="\t", file=f)
            print("PEMap_rate", oj.PEMap_rate, sep="\t", file=f)
            print("insert_rate", oj.insert_rate, sep="\t", file=f)

        

if __name__ == "__main__":
    main()

#bed_file = "/project/BIT/keke.zou/homebin/cbNIPT_con/PH_SCv1.bed"
#bam_file = "/project/NIPT/G161002DNGS_CBNIPT/170905_Paternity/20180227/PB2017082901M/PB2017082901M.bam"
#bam_file = "/project/PGS_PGD/Z1707ACNGS_UniSeqPGS/20180831/cut_len/CJI1112_S24/CJI1112_S24.bam"
# bam_file = "/project/BIT/keke.zou/analysis/chongqingfuyou/20180929/R18041624LU01-YJ/R18041624LU01-YJ.bam"
# bed_file = "/data/Customer/chongqifuyou/bed/S07604514_Covered.bed"
# oj = BamQc(bam_file,bed_file)
# oj.get_stat()
#


