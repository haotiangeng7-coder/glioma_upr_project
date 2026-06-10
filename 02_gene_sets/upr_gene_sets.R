###############################################################################
# upr_gene_sets.R
# 定义UPR三条信号臂及广谱基因集
# 来源：MSigDB, GO, KEGG, 文献整理
###############################################################################

source("00_setup/config.R")
library(msigdbr)

# =============================================================================
# 1. 手动整理的UPR核心基因（文献来源）
# =============================================================================

# --- IRE1α-XBP1通路 ---
# IRE1α(ERN1)激活后剪切XBP1 mRNA产生XBP1s，驱动ER伴侣蛋白和ERAD基因表达
IRE1_XBP1_genes <- c(
  # 核心信号分子
  "ERN1",       # IRE1α，ER跨膜激酶/核酸内切酶

"XBP1",       # XBP1，剪切后激活转录
  "DNAJB9",     # ERdj4，XBP1s靶基因，ER伴侣蛋白
  "EDEM1",      # ERAD相关，降解错误折叠蛋白
  "HERPUD1",    # ERAD组分
  "SEC61A1",    # 蛋白转运通道
  "SEC61B",     # 蛋白转运通道亚基
  "BLOC1S1",    # XBP1s调控靶基因
  "SYVN1",      # HRD1, E3泛素连接酶，ERAD
  "DERL1",      # Derlin-1, ERAD
  "DERL2",      # Derlin-2, ERAD
  "OS9",        # ERAD lectin
  "SEL1L",      # ERAD cofactor
  "VCP",        # p97, ERAD AAA ATPase
  "PDIA6",      # 蛋白二硫键异构酶
  "HYOU1",      # GRP170/ORP150, ER伴侣蛋白
  "DNAJC3",     # p58IPK, 负调控PERK
  "SERP1",      # RAMP4, 应激相关膜蛋白
  "KDELR1",     # ER retention receptor
  "SSR1"        # TRAP alpha
)

# --- PERK-eIF2α-ATF4通路 ---
# PERK(EIF2AK3)磷酸化eIF2α，全局翻译抑制但选择性翻译ATF4
PERK_ATF4_genes <- c(
  # 核心信号分子
  "EIF2AK3",    # PERK，ER跨膜激酶
  "EIF2S1",     # eIF2α，翻译起始因子
  "ATF4",       # 转录因子，ISR核心效应子
  "DDIT3",      # CHOP/GADD153，促凋亡转录因子
  "ASNS",       # 天冬酰胺合成酶，ATF4靶基因
  "TRIB3",      # ATF4/CHOP靶基因，抑制Akt
  "PPP1R15A",   # GADD34，去磷酸化eIF2α（负反馈）
  "ATF3",       # 应激诱导转录因子
  "CEBPB",      # C/EBP beta
  "SLC7A11",    # xCT, 胱氨酸/谷氨酸转运体
  "PHGDH",      # 丝氨酸合成
  "PSAT1",      # 丝氨酸合成
  "SLC7A5",     # 氨基酸转运体
  "SLC7A1",     # 氨基酸转运体
  "VEGFA",      # ATF4靶基因，血管生成
  "CEBPG",      # C/EBP gamma
  "HERPUD1",    # 也受PERK调控
  "SESN2",      # Sestrin 2, 应激响应
  "WARS1",      # tRNA synthetase, ISR靶基因
  "CHAC1"       # 谷胱甘肽降解, CHOP靶基因
)

# --- ATF6通路 ---
# ATF6在ER应激下转运至高尔基体被S1P/S2P切割激活
ATF6_genes <- c(
  # 核心信号分子
  "ATF6",       # ATF6α，ER跨膜转录因子
  "ATF6B",      # ATF6β
  "HSPA5",      # GRP78/BiP, ER主伴侣蛋白, 调控UPR三个传感器
  "HSP90B1",    # GRP94, ER伴侣蛋白
  "PDIA4",      # 蛋白二硫键异构酶
  "CALR",       # 钙网蛋白，蛋白质折叠质控
  "CANX",       # 连接蛋白，蛋白质折叠质控
  "PDIA3",      # ERp57
  "P4HB",       # PDI, 蛋白二硫键异构酶
  "MBTPS1",     # S1P, ATF6切割酶
  "MBTPS2",     # S2P, ATF6切割酶
  "CREB3",      # OASIS, 类ATF6转录因子
  "CREB3L1",    # 类ATF6
  "CREB3L2",    # 类ATF6
  "SDF2L1",     # ER质控
  "PPIB",       # Cyclophilin B
  "FKBP2",      # ER FKBP
  "DNAJB11",    # ERdj3, ER伴侣蛋白
  "CRELD2",     # ATF6靶基因
  "MANF"        # 中脑星形胶质细胞衍生神经营养因子
)

# --- 综合UPR广谱基因集 ---
# 合并三条通路 + 额外UPR相关基因
UPR_additional <- c(
  # ER应激传感器共同调控
  "ERO1A",      # ER氧化酶
  "ERO1B",      # ER氧化酶
  "UGGT1",      # UDP-glucose glycoprotein glucosyltransferase
  "UGGT2",
  "RPN1",       # 核糖体蛋白，N-糖基化
  "RPN2",
  "STT3A",      # 寡糖基转移酶
  "STT3B",
  "DDOST",      # 寡糖基转移酶
  "DAD1",       # 寡糖基转移酶
  "DNAJB1",     # HSP40家族
  "DNAJB2",
  "HSPA1A",     # HSP70
  "HSPA1B",
  "HSPA8",      # HSC70
  "BAX",        # 促凋亡（UPR持续激活的下游）
  "BCL2",       # 抗凋亡
  "CASP4",      # ER应激特异性caspase
  "CASP12"      # ER应激特异性caspase（人类中假基因化）
)

UPR_broad_genes <- unique(c(
  IRE1_XBP1_genes,
  PERK_ATF4_genes,
  ATF6_genes,
  UPR_additional
))

# =============================================================================
# 2. 从MSigDB获取GO UPR基因集（补充验证）
# =============================================================================

message("Fetching UPR gene sets from MSigDB...")

# GO:BP - 所有UPR相关GO terms
msigdb_go <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
upr_go_terms <- msigdb_go[grep("UNFOLDED_PROTEIN|PERK_MEDIATED|IRE1_MEDIATED|ATF6_MEDIATED",
                                msigdb_go$gs_name), ]
# 合并所有UPR相关GO基因
upr_go <- upr_go_terms

# Hallmark UPR
msigdb_h <- msigdbr(species = "Homo sapiens", collection = "H")
upr_hallmark <- msigdb_h[msigdb_h$gs_name == "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", ]

# KEGG protein processing in ER (may not be available in KEGG_LEGACY)
tryCatch({
  msigdb_kegg <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
  er_kegg <- msigdb_kegg[grep("ENDOPLASMIC|PROTEIN_EXPORT", msigdb_kegg$gs_name), ]
}, error = function(e) {
  er_kegg <<- data.frame(gene_symbol = character(0))
})

# Reactome UPR
msigdb_reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
upr_reactome <- msigdb_reactome[grep("UNFOLDED_PROTEIN", msigdb_reactome$gs_name), ]

# =============================================================================
# 3. 整合所有基因集
# =============================================================================

upr_gene_list <- list(
  IRE1_XBP1 = IRE1_XBP1_genes,
  PERK_ATF4 = PERK_ATF4_genes,
  ATF6      = ATF6_genes,
  UPR_broad = UPR_broad_genes,
  GO_UPR    = unique(upr_go$gene_symbol),
  Hallmark_UPR = unique(upr_hallmark$gene_symbol),
  KEGG_ER_processing = unique(er_kegg$gene_symbol),
  Reactome_UPR = unique(upr_reactome$gene_symbol)
)

# 打印基因集大小
message("\n=== UPR Gene Set Summary ===")
for (nm in names(upr_gene_list)) {
  message(sprintf("  %-25s: %d genes", nm, length(upr_gene_list[[nm]])))
}

# =============================================================================
# 4. 免疫相关基因集（用于后续分析）
# =============================================================================

immune_checkpoint_genes <- c(
  "PDCD1",    # PD-1
  "CD274",    # PD-L1
  "PDCD1LG2", # PD-L2
  "CTLA4",    # CTLA-4
  "HAVCR2",   # TIM-3
  "LAG3",     # LAG-3
  "TIGIT",    # TIGIT
  "VSIR",     # VISTA
  "IDO1",     # IDO1
  "CD276",    # B7-H3
  "VTCN1",    # B7-H4
  "TNFRSF9",  # 4-1BB
  "ICOS",     # ICOS
  "TNFRSF4",  # OX40
  "TNFRSF18", # GITR
  "CD40",     # CD40
  "CD40LG",   # CD40L
  "CD80",     # B7-1
  "CD86",     # B7-2
  "TNFSF4",   # OX40L
  "TNFSF9",   # 4-1BBL
  "ADORA2A",  # A2AR
  "BTLA",     # BTLA
  "SIGLEC15"  # Siglec-15
)

# T细胞耗竭标志基因
t_exhaustion_genes <- c(
  "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4",
  "ENTPD1", "LAYN", "TOX", "TOX2", "CXCL13",
  "BATF", "IRF4", "NFATC1", "NR4A1", "NR4A2"
)

# M1/M2巨噬细胞标志
m1_markers <- c("TNF", "IL1B", "IL6", "NOS2", "CD80", "CD86", "IRF5", "STAT1")
m2_markers <- c("ARG1", "MRC1", "CD163", "IL10", "TGFB1", "CCL22", "IRF4", "STAT6")

# =============================================================================
# 5. 保存基因集
# =============================================================================

save(
  upr_gene_list,
  IRE1_XBP1_genes,
  PERK_ATF4_genes,
  ATF6_genes,
  UPR_broad_genes,
  immune_checkpoint_genes,
  t_exhaustion_genes,
  m1_markers,
  m2_markers,
  file = file.path(DATA_PROC, "upr_gene_sets.RData")
)

message("\nGene sets saved to: ", file.path(DATA_PROC, "upr_gene_sets.RData"))
