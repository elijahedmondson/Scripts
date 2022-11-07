setwd("/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/Background")
load(file = "/Users/elijahedmondson/Desktop/R/QTL/WD/GRSD_master.Rdata")
pheno = data.frame(row.names = Background$row.names, 
                   sex = as.numeric(Background$sex == "M"),  
                   LSA = as.numeric(Background$Lymphoma),
                   NN = as.numeric(Background$Non.neoplastic),
                   PulACA = as.numeric(Background$Pulmonary.Adenocarcinoma),
                   PulAd = as.numeric(Background$Pulmonary.Adenoma),
                   PulSrc = as.numeric(Background$Pulmonary.Sarcomatoid.Carcinoma),
                   BHGT = as.numeric(Background$Bilateral.Harderian.Gland.Tumors),
                   CPT = as.numeric(Background$Choroid.Plexus.Tumor),
                   Dpolyp = as.numeric(Background$Duodenal.Polyp),
                   EndoPolyp = as.numeric(Background$Endometrial.Stromal.Polyp),
                   Ectoderm = as.numeric(Background$Ectoderm),
                   Endoderm = as.numeric(Background$Endoderm),
                   Ependymoma = as.numeric(Background$Ependymoma),
                   FSA = as.numeric(Background$Fibrosarcoma),
                   GastSCC = as.numeric(Background$Gastric.Squamous.Cell.Carcinoma),
                   GCT = as.numeric(Background$Granulosa.Cell.Tumor),
                   HGAca = as.numeric(Background$Harderian.Gland.Adenocarcinoma),
                   HGAd = as.numeric(Background$Harderian.Gland.Adenoma),
                   Hemangioma = as.numeric(Background$Hemangioma),
                   HSA = as.numeric(Background$Hemangiosarcoma),
                   Hepatoblastoma = as.numeric(Background$Hepatoblastoma),
                   HCC = as.numeric(Background$Hepatocellular.Carcinoma),
                   HCAd = as.numeric(Background$Hepatocellular.Adenoma),
                   HS = as.numeric(Background$Histiocytic.Sarcoma),
                   IntestinalACA = as.numeric(Background$Intestinal.Adenocarcinoma),
                   IntestinalNeo = as.numeric(Background$Intestinal.Neoplasms),
                   Intracranial = as.numeric(Background$Intracranial.Tumors),
                   Islet = as.numeric(Background$Islet.Cell.Carcinoma),
                   Leiomyosarcoma = as.numeric(Background$Leiomyosarcoma),
                   MalMammary = as.numeric(Background$Malignant.Mammary.Tumors),
                   MalOvarian = as.numeric(Background$Malignant.Ovarian.Tumors),
                   MammAdenoacanthoma = as.numeric(Background$Mammary.Gland.Adenoacanthoma),
                   MammAdenocarcinoma = as.numeric(Background$Mammary.Gland.Adenocarcinoma),
                   Meningioma = as.numeric(Background$Meningioma),
                   Mesoderm = as.numeric(Background$Mesoderm),
                   Mesothelioma = as.numeric(Background$Mesothelioma),
                   Met = as.numeric(Background$metastasis),
                   Metastatic = as.numeric(Background$Metastatic.Tumors),
                   Myeloid.Leukemia = as.numeric(Background$Myeloid.Leukemia),
                   Myxosarcoma = as.numeric(Background$Myxosarcoma),
                   Odontogenic = as.numeric(Background$Odontogenic.Tumor),
                   Osteoma = as.numeric(Background$Osteoma),
                   OSA = as.numeric(Background$Osteosarcoma),
                   OvarianCarc = as.numeric(Background$Ovarian.Carcinoma),
                   Pheo = as.numeric(Background$Pheochromocytoma),
                   PitAd = as.numeric(Background$Pituitary.Adenoma),
                   PulAd = as.numeric(Background$Pulmonary.Adenoma),
                   PulMetastasis = as.numeric(Background$Pulmonary.Metastases),
                   RCC = as.numeric(Background$Renal.Cell.Carcinoma),
                   RhSA = as.numeric(Background$Rhabdomyosarcoma),
                   STS = as.numeric(Background$Soft.Tissue.Sarcomas),
                   ThyroidAd = as.numeric(Background$Thyroid.Adenoma),
                   ThyroidCarc = as.numeric(Background$Thyroid.Carcinoma),
                   TubulostromalACA = as.numeric(Background$Tubulostromal.Adenocarcinoma),
                   TubulostromalAd = as.numeric(Background$Tubulostromal.Adenoma),
                   UndiffSarc = as.numeric(Background$Undifferentiated.Sarcoma),
                   UterineStromalSarc = as.numeric(Background$Uterine.Stromal.Sarcoma),
                   UnilHarderian = as.numeric(Background$Unilateral.Harderian.Gland.Tumors),
                   Epidermal = as.numeric(Background$Tumors.of.Epidermis))
covar = data.frame(sex = as.numeric(Background$sex == "M"))
addcovar = covar
rownames(covar) = rownames(pheno)
rownames(addcovar) = rownames(pheno)

LM.qtl = scanone(pheno = pheno, pheno.col = "LSA", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.LSA.Rdata")
plot(LM.qtl)
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "LSA", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.LSA.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "NN", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.NN.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "NN", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.NN.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "PulACA", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.PulmonaryAdenocarcinoma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "PulACA", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.PulmonaryAdenocarcinoma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "PulAd", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.PulAdenoma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "PulAd", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.PulAdenoma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "PulSrc", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.PulSarcomatoidCarc.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "PulSrc", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.PulSarcomatoidCarc.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "BHGT", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.BHGT.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "BHGT", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.BHGT.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "CPT", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.CPT.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "CPT", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.CPT.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Dpolyp", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Dpolyp,Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Dpolyp", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Dpolyp.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "EndoPolyp", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.EndoPolyp.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "EndoPolyp", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.EndoPolyp.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Ependymoma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Ependymoma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Ependymoma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Ependymoma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "FSA", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.FSA.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "FSA", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.FSA.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "GastSCC", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.GastSCC.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "GastSCC", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.GastSCC.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "GCT", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.GCT.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "GCT", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.GCT.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "HGAca", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.HGAca.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "HGAca", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.HGAca.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "HGAd", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.HGAd.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "HGAd", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.HGAd.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Hemangioma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Hemangioma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Hemangioma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Hemangioma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "HSA", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.HSA.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "HSA", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.HSA.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Hepatoblastoma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Hepatoblastoma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Hepatoblastoma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Hepatoblastoma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "HCC", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.HCC.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "HCC", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.HCC.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "HCAd", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.HCAd.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "HCAd", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.HCAd.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "HS", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.HS.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "HS", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.HS.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "IntestinalACA", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.IntestinalACA.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "IntestinalACA", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.IntestinalACA.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "IntestinalNeo", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.IntestinalNeo.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "IntestinalNeo", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.IntestinalNeo.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Intracranial", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Intracranial.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Intracranial", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Intracranial.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Islet", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Islet.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Islet", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Islet.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Myeloid.Leukemia", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Myeloid.Leukemia.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Myeloid.Leukemia", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Myeloid.Leukemia.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Myxosarcoma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Myxosarcoma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Myxosarcoma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Myxosarcoma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Odontogenic", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Odontogenic.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Odontogenic", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Odontogenic.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Osteoma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Osteoma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Osteoma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Osteoma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "OSA", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.OSA.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "OSA", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.OSA.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Leiomyosarcoma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Leiomyosarcoma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Leiomyosarcoma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Leiomyosarcoma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "MalMammary", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.MalMammary.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "MalMammary", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.MalMammary.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "MalOvarian", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.MalOvarian.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "MalOvarian", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.MalOvarian.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "MammAdenoacanthoma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.MammAdenoacanthoma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "MammAdenoacanthoma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.MammAdenoacanthoma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "MammAdenocarcinoma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.MammAdenocarcinoma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "MammAdenocarcinoma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.MammAdenocarcinoma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Meningioma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Meningioma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Meningioma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Meningioma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Mesothelioma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Mesothelioma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Mesothelioma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Mesothelioma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Met", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Met.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Met", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Met.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Metastatic", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Metastatic.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Metastatic", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Metastatic.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Ectoderm", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Ectoderm.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Ectoderm", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Ectoderm.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Endoderm", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Endoderm.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Endoderm", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Endoderm.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Mesoderm", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Mesoderm.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Mesoderm", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Mesoderm.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "OvarianCarc", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.OvarianCarc.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "OvarianCarc", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.OvarianCarc.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Pheo", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Pheo.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Pheo", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Pheo.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "PitAd", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.PitAd.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "PitAd", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.PitAd.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "PulAd", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.PulAd.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "PulAd", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.PulAd.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "PulMetastasis", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.PulMetastasis.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "PulMetastasis", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.PulMetastasis.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "RCC", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.RCC.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "RCC", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.RCC.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "RhSA", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.RhSA.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "RhSA", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.RhSA.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "STS", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.STS.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "STS", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.STS.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "ThyroidAd", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.ThyroidAd.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "ThyroidAd", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.ThyroidAd.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "ThyroidCarc", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.ThyroidCarc.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "ThyroidCarc", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.ThyroidCarc.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "TubulostromalACA", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.TubulostromalACA.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "TubulostromalACA", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.TubulostromalACA.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "TubulostromalAd", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.TubulostromalAd.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "TubulostromalAd", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.TubulostromalAd.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "UndiffSarc", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.UndiffSarc.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "UndiffSarc", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.UndiffSarc.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "UterineStromalSarc", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.UterineStromalSarc.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "UterineStromalSarc", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.UterineStromalSarc.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "UnilHarderian", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.UnilHarderian.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "UnilHarderian", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.UnilHarderian.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Epidermal", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Epidermal.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Epidermal", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Epidermal.Rdata")
rm(AM.qtl, LM.qtl)