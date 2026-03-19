library(ggplot2)
library(extrafont)
library(systemfonts)
library(showtext)
library(ggtext)
library(patchwork)
library(ggfx)
library(ggforce)

## 1. Typography Setup ----
# Load Michroma from Google Fonts
font_add_google("Space Grotesk", "space")
showtext_auto()
showtext::showtext_opts(dpi = 600)

# These commands help R find and manage the system fonts
font_paths()  
font_files()
font_families()

# Automatically triggers showtext to ensure custom fonts 
# render correctly in high-res PNG outputs
trace(grDevices::png, exit = quote({
  showtext::showtext_begin()
}), print = FALSE)

kraken_palette <- c(
  "Macro Complexity"          = "#000E3B",
  "Trophic Integration"       = "#68A2B9",
  "Energy Transport"          = "#E9072B",
  "Trophic Role Asymmetry"    = "#99D9D9",
  "Control Heterogeneity"     = "#005A5B",
  "Basal Control"             = "#9BBF80"
)

# Converts the named vector into a dataframe for easier mapping in ggplot
pal_df <- tibble(
  label   = names(kraken_palette), # printed label
  colour  = kraken_palette, # colour code
  value   = factor(1:length(kraken_palette)) # module number
)


# continuous ramp
seattle_anchors <- c("#C0DADC", "#93C1D2", "#638596", "#324B5C", "#000E3B")

# 2. Create the Generator Function
# This creates a function that can interpolate any number of colors
seattle_abyssal_gen <- colorRampPalette(seattle_anchors)

# diverging colour
pal_diverge <- tibble(
  low  = "#005A9C",  # darker trench
  mid  = "#F5F5F5",
  high = "#00843D"   # brighter copper
)

seattle_div <- c("#001457", "#638596", "#F0F0F0", "#86BD26", "#FFB612")
seattle_div <- c("#E9072B", "#F5F5F5", "#68A2B9")

pal_diverge_gen <- colorRampPalette(seattle_div)

# additional categorical colours

secondary_palette <- c(
  "#AFD2E3", 
  "#283339", 
  "#355464"
  )

# stability palette
stability_palette <- c(
  "#198566", 
  "#771985", 
  "#856819", 
  "#FFB612"
)

pca_col <- secondary_palette[3]

figure_theme <- function() {
  theme_bw() %+replace% 
    theme(
      # Text & Axes
      text = element_text(family = "space", color = "#001628"),
      plot.title = element_text(face = "bold", size = rel(1.3), margin = margin(b = 10)),
      axis.text = element_text(size = rel(0.7), color = "#001628"),
      axis.title = element_text(size = rel(0.9), face = "bold"),
      
      # Grid & Background
      panel.grid.major = element_line(color = "#A5ACAF", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      
      # Border & Strip Alignment
      panel.border = element_rect(color = "#001628", fill = NA, linewidth = 1),
      strip.background = element_rect(color = "#001628", fill = "white", linewidth = 1),
      
      # No gap between strip and panel
      panel.spacing = unit(0, "pt"),
      
      # Legend
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
}
