using FlowAtlas
channelMap = Dict([

    "FJComp-355 379_28-A" => "CD3-IgD", 
    "CD3" => "CD3-IgD"

    "FJComp-355 560_40-A" => "CD8", 
    "FJComp-405 585_15-A" => "CD19",

    # "FJComp-355 820_60-A" => "CD4",
    "FJComp-355 670_30-A" => "CD4",

    "FJComp-640 780_60-A" => "CCR7",
    "FJComp-405 780_60-A" => "CD45RA", 

    "FJComp-561 780_60-A" => "CD127", 
    "FJComp-640 670_30-A" => "CD25", 

    "FJComp-561 610_20-A" => "Helios", 
    "FJComp-561 585_15-A" => "Foxp3-IgM", 
    "Foxp3" => "Foxp3-IgM",

    "FJComp-405 710_40-A" => "PD-1", 
    "FJComp-640 730_35-A" => "CXCR5", 

    "FJComp-405 670_30-A" => "CCR6", 
    "FJComp-488 715_30-A" => "CXCR3", 

    "FJComp-405 605_40-A" => "CCR4", 
    "FJComp-488 525_50-A" => "CCR10", 

    "FJComp-405 450_50-A" => "CD103", 
    "FJComp-355 740_35-A" => "CD69",
    "FJComp-405 515_20-A" => "HLA-DR"
])

FlowAtlas.run("data/workspace.wsp",glob"data/*/*.fcs",channelMap=channelMap)