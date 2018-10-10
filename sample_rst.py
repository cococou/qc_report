#!/usr/bin/env python3
__author__ = 'zoukeke'

import os
import shutil
import argparse

def get_html(stat_dict=None):#,sample=None,atcg_pic=None,QC_pic=None,insize_pic=None,heatmap=None):

    sample_html = """
{sample}
{main_line}

QC report of bam
----------------

======================  =======================
参数                     值
======================  =======================
total bases             {bases}
total reads             {reads}
reads length            {len_reads}
insert size median      {insert_mid}
Q20 rate                {Q20_rate}
Q30 rate                {Q30_rate}
GC rate                 {gc_rate}
GC align rate           {gc_align_rate}
GC unique rate          {gc_uniq_rate}
mismatch rate           {mismatch_rate}
deletion rate           {deletion_rate}
insert rate             {insert_rate}
gap rate                {gap_rate}
mapping quality mean    {map_quali}
mapping rate            {mapping_rate}
unique mapping rate     {UniqueMap_rate}
duplication rate        {duplication_rate}
PEMap rate              {PEMap_rate}
UniqueMap rate          {UniqueMap_rate}
======================  =======================


QC report of targeted bam
-------------------------

==================================  =======================
参数                                 值
==================================  =======================
target length                       {target_len}
coverage rate                       {cover_rate}
target region base's mean depth     {target_mean_depth}
more than 20% mean depth rate       {cover_cond_rate_20}
more than 50% mean depth rate       {cover_cond_rate_50}
==================================  =======================

位置QC质量分布图
-----------------
.. image:: _static/{QC_pic}

位置ATCG含量分布图
-------------------
.. image:: _static/{atcg_pic}

insert size图
---------------
.. image:: _static/{insize_pic}

region深度热度图
------------------
.. image:: _static/{heatmap}

""".format(bases=stat_dict['bases'],reads=stat_dict['reads'],len_reads=stat_dict['len_reads'],insert_mid=stat_dict['insert_mid'],Q20_rate=stat_dict['Q20_rate'],Q30_rate=stat_dict['Q30_rate'],gc_rate=stat_dict['gc_rate'],gc_align_rate=stat_dict['gc_align_rate'],gc_uniq_rate=stat_dict['gc_uniq_rate'],mismatch_rate=stat_dict['mismatch_rate'],deletion_rate=stat_dict['deletion_rate'],insert_rate=stat_dict['insert_rate'],map_quali=stat_dict['map_quali'],mapping_rate=stat_dict['mapping_rate'],UniqueMap_rate=stat_dict['UniqueMap_rate'],duplication_rate=stat_dict['duplication_rate'],PEMap_rate=stat_dict['PEMap_rate'],target_len=stat_dict['target_len'],cover_rate=stat_dict['cover_rate'],target_mean_depth=stat_dict['target_mean_depth'],cover_cond_rate_20=stat_dict['cover_cond_rate_20'],cover_cond_rate_50=stat_dict['cover_cond_rate_50'],sample=stat_dict['sample'],main_line=len(stat_dict['sample'],)*"=",atcg_pic=stat_dict['atcg_pic'],QC_pic=stat_dict['QC_pic'],insize_pic=stat_dict['insize_pic'],heatmap=stat_dict['heatmap'],gap_rate=stat_dict['gap_rate'])
    return sample_html

def get_info(f=None):
    stat_dict = {'bases': 'NA', 'reads': 'NA', "len_reads": "NA", "insert_mid": "NA", 'Q20_rate': "NA",
                 "Q30_rate": "NA", "gc_rate": "NA", 'gc_align_rate': "NA", "gc_uniq_rate": "NA", 'mismatch_rate': "NA",
                 'deletion_rate': "NA", 'insert_rate': "NA", 'gap_rate': "NA", 'map_quali': "NA", 'mapping_rate': "NA",
                 'UniqueMap_rate': "NA", 'duplication_rate': "NA", 'PEMap_rate': "NA", 'UniqueMap_rate': "NA",
                 'target_len': "NA", 'cover_rate': "NA", 'target_mean_depth': "NA", 'cover_cond_rate_20': "NA",
                 'cover_cond_rate_50': "NA",'sample':'NA','heatmap':'NA',"atcg_pic":"NA","QC_pic":'NA','insize_pic':"NA"}
    with open(f) as ff:
        for line in ff:
            line = line.strip()
            if not line:continue
            key,value = line.split("\t")[0:2]
            if key in stat_dict:
                if "_rate" in key:
                    stat_dict[key] = "%.2f%%" % (float(value) * 100)
                elif "mean" in key:
                    stat_dict[key] = "%.2f" % float(value)
                else:
                    stat_dict[key] = value
    return stat_dict

def copy_png(qc_dir=None,pic_dir=None):
    fs = [os.path.join(qc_dir,i) for i in os.listdir(qc_dir) if i.endswith(".png")]
    for x in fs:
        if os.path.isfile(x):
            shutil.copyfile(x,os.path.join(pic_dir,os.path.basename(x))) #覆盖拷贝

def get_sample_dir(rawdir=None):
    fs = os.listdir(rawdir)
    qc_dir = []
    for f in fs:
        ff = os.path.join(rawdir,f,"qc")
        if os.path.exists(ff):
            ff_fs = os.listdir(ff)
            if "stat.txt" in ff_fs:
                qc_dir.append(ff)
    return qc_dir

def regenerate_temp(outdir=None,template=None):
    os.makedirs(outdir,exist_ok=True)
    os.chdir(outdir)
    try:
        shutil.copyfile(template,os.path.join(outdir,os.path.basename(template)))
        os.system("tar -zxvf {}".format(os.path.join(outdir,os.path.basename(template))))
    except:
        os.system("git clone {}".format(template))

def get_static_pic(qc_dir=None,stat_dict=None):
    pics = os.listdir(qc_dir)
    sample = os.path.basename(os.path.dirname(qc_dir)).split('.')[0]
    stat_dict['sample'] = sample
    #,sample=None,atcg_pic=None,QC_pic=None,insize_pic=None,heatmap=None))
    #PB2017082901M_gc.png  PB2017082901M_insert.png  PB2017082901M_quality.png  PB2017082901M_target.png
    for pic in pics:
        if pic.endswith("_gc.png"):
            stat_dict['atcg_pic'] = pic
        elif pic.endswith("_quality.png"):
            stat_dict['QC_pic'] = pic
        elif pic.endswith("_insert.png"):
            stat_dict["insize_pic"] = pic
        elif pic.endswith("_target.png"):
            stat_dict['heatmap'] = pic
    return stat_dict

def ch_index_rst(index_rst_file=None,sample_rsts=None):
    index_rst = """
.. report documentation master file, created by
   sphinx-quickstart on Wed Oct 10 10:55:28 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to QC reports!
======================

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   {samples}

search
==================

* :ref:`search`

    """.format(samples=sample_rsts)

    with open(index_rst_file,'w') as f:
        print(index_rst,file=f)


def main(outdir=None,rawdir=None,report_dir=None,template=None):
    qc_dirs =get_sample_dir(rawdir) #得到有结果的qc文件的路径 [sample_dir1,sample_dir2]
    regenerate_temp(outdir, template)  #生成模板
    pic_dir=os.path.join(outdir,"_static") #存放图片的路径
    if not report_dir:
        report_dir=os.path.join(outdir,'report')

    #x= list(map(lambda qc_dir:copy_png(qc_dir, pic_dir),qc_dirs)) #移动文件
    sample_rsts = []
    for qc_dir in qc_dirs:
        copy_png(qc_dir, pic_dir) #准备图片文件
        stat_dict = get_info(os.path.join(qc_dir,"stat.txt"))
        stat_dict = get_static_pic(qc_dir, stat_dict) # get png
        html = get_html(stat_dict)
        sample_rst = os.path.join(outdir,stat_dict['sample']+".rst")
        with open(sample_rst,'w') as f:
            print(html,file=f)
            sample_rsts.append(sample_rst)
    sample_rsts = '\n\n   '.join([os.path.basename(i) for i in sample_rsts])
    index_rst_file = os.path.join(outdir,"index.rst")
    ch_index_rst(index_rst_file,sample_rsts) #生成index.rst 文件
    os.chdir(outdir)
    os.system("sphinx-build {outdir} {report_dir}".format(outdir=outdir,report_dir=report_dir))


def get_par():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r","--rawdir",help="samples' directory",default=None)
    parser.add_argument("-o","--outdir",help="template directory, default .",default=os.getcwd())
    parser.add_argument("-t","--template",help="raw template file or git file,default is /project/BIT/keke.zou/pycharm/p3/sphinx/report_v1.tar.gz",default="/project/BIT/keke.zou/pycharm/p3/sphinx_html/report_v1.tar.gz")
    parser.add_argument("-p","--report_dir",help="report directory,default is current report",default=os.path.join(os.getcwd(),"report"))
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    par = get_par()
    main(outdir=os.path.abspath(par.outdir),rawdir=os.path.abspath(par.rawdir),template=par.template,report_dir=os.path.abspath(par.report_dir))
# outdir="/project/BIT/keke.zou/pycharm/p3/sphinx/test"
# rawdir = "/project/BIT/keke.zou/pycharm/p3/"
# main(outdir,rawdir)


