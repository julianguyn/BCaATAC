hallmark_categories <- data.frame(
  pathway = c(
    "ADIPOGENESIS","ALLOGRAFT_REJECTION","ANDROGEN_RESPONSE","ANGIOGENESIS",
    "APICAL_JUNCTION","APICAL_SURFACE","APOPTOSIS","BILE_ACID_METABOLISM",
    "CHOLESTEROL_HOMEOSTASIS","COAGULATION","COMPLEMENT","DNA_REPAIR",
    "E2F_TARGETS","EPITHELIAL_MESENCHYMAL_TRANSITION","ESTROGEN_RESPONSE_EARLY",
    "ESTROGEN_RESPONSE_LATE","FATTY_ACID_METABOLISM","G2M_CHECKPOINT",
    "GLYCOLYSIS","HEDGEHOG_SIGNALING","HEME_METABOLISM","HYPOXIA",
    "IL2_STAT5_SIGNALING","IL6_JAK_STAT3_SIGNALING","INFLAMMATORY_RESPONSE",
    "INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE","KRAS_SIGNALING_DN",
    "KRAS_SIGNALING_UP","MITOTIC_SPINDLE","MTORC1_SIGNALING","MYC_TARGETS_V1",
    "MYC_TARGETS_V2","MYOGENESIS","NOTCH_SIGNALING","OXIDATIVE_PHOSPHORYLATION",
    "P53_PATHWAY","PANCREAS_BETA_CELLS","PEROXISOME","PI3K_AKT_MTOR_SIGNALING",
    "PROTEIN_SECRETION","REACTIVE_OXYGEN_SPECIES_PATHWAY","SPERMATOGENESIS",
    "TGF_BETA_SIGNALING","TNFA_SIGNALING_VIA_NFKB","UNFOLDED_PROTEIN_RESPONSE",
    "UV_RESPONSE_DN","UV_RESPONSE_UP","WNT_BETA_CATENIN_SIGNALING",
    "XENOBIOTIC_METABOLISM"
  ),
  category = c(
    "Development","Immune","Hormone_Response","Development",
    "Development","Development","Cellular_Stress","Metabolism",
    "Metabolism","Immune","Immune","DNA_Damage",
    "Proliferation","Development","Hormone_Response",
    "Hormone_Response","Metabolism","Proliferation",
    "Metabolism","Signaling","Metabolism","Cellular_Stress",
    "Immune","Immune","Immune",
    "Immune","Immune","Signaling",
    "Signaling","Proliferation","Signaling","Proliferation",
    "Proliferation","Development","Signaling","Metabolism",
    "Proliferation","Development","Metabolism","Signaling",
    "Development","Cellular_Stress","Development",
    "Signaling","Immune","Cellular_Stress",
    "Cellular_Stress","Cellular_Stress","Signaling",
    "Metabolism"
  ),
  stringsAsFactors = FALSE
)