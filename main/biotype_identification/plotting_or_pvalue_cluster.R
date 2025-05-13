abc = calculate_stat_for_res(collect_res,num_it)

scz = "PGC010"

cognitive_impairment = c("GCST90309368",
                         "GCST90179115",
                         "GCST90162549",
                         "GCST90162548",
                         "GCST90162547",
                         "GCST90179117")


neuroticism = c("BIG01",
                "GCST006478",
                "GCST006948",
                "GCST006952",
                "GCST90029028",
                "GCST90041871",
                "GCST90041881",
                "GCST90042775",
                "GCST90271701",
                "GCST90271707",
                "GCST90271708",
                "GCST90041879",
                "GCST90042756",
                "GCST90271710")

depressed = c("GCST005902",
              "GCST009982",
              "GCST90014431",
              "GCST90014435",
              "GCST90042682",
              "GCST90042772",
              "GCST90042804",
              "GCST90042861",
              "GCST90042862",
              "GCST90044036",
              "GCST90044037",
              "GCST90267276",
              "GCST90309343",
              "PGC013",
              "PGC032",
              "GCST90271709",
              "GCST90042802",
              "GCST90014436",
              "GCST90042695",
              "PGC033")

DISCOVERY_DATA = DiscoveryRES     # DELETE
VALIDATION_DATA = ValidationRES     # DELETE

best_prs = fread("/data/projects/on_going/PGScope/resources/best_prs.tsv", data.table = T)

discovery = "#DDEAF9"
validation ="#F9F2DD"

############################TRAINING###################################################

mean_plot(df_discovery, DISCOVERY_DATA,features_sizes, "", best_prs, scz, plot_b =discovery,y_min = -0.5,y_max = 1.5) 
mean_plot(df_discovery, DISCOVERY_DATA,features_sizes,"", best_prs, cognitive_impairment, plot_b =discovery,y_min =-0.6,y_max =0.25)
mean_plot(df_discovery, DISCOVERY_DATA,features_sizes, "", best_prs, neuroticism, plot_b =discovery,y_min = -0.8,y_max =1.3) 
mean_plot(df_discovery,DISCOVERY_DATA,features_sizes,"", best_prs, depressed, plot_b =discovery, y_min =-0.5,y_max =1) 

or_plot("Discovery",features_sizes, plot_b =discovery) 
pvalue_plot("Discovery",features_sizes, plot_b =discovery) 

##############################TEST-SCZ####################################################

mean_plot(df_validation,VALIDATION_DATA,features_sizes, "", best_prs, scz, plot_b=validation,y_min =-0.5,y_max =1.5)
mean_plot(df_validation,VALIDATION_DATA,features_sizes,"", best_prs, cognitive_impairment, plot_b=validation,y_min =-0.6,y_max =0.25) 
mean_plot(df_validation,VALIDATION_DATA,features_sizes,"", best_prs, neuroticism, plot_b=validation, y_min =-0.8,y_max =1.3) 
mean_plot(df_validation,VALIDATION_DATA,features_sizes,"", best_prs, depressed, plot_b=validation, y_min =-0.5,y_max =1)

or_plot("Validation",features_sizes, plot_b=validation)
pvalue_plot("Validation",features_sizes, plot_b=validation) 

