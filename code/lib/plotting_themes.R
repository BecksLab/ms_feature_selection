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
  "Macro Complexity"       = "#001628", # Deep Sea Blue (Kraken)
  "Trophic Integration"    = "#69BE28", # Action Green (Seahawks)
  "Energy Transport"       = "#FBD800", # Storm Gold (Storm)
  "Trophic Asymmetry"      = "#005A5B", # Torrent Teal (Torrent)
  "Control Heterogeneity"  = "#E9072B", # Red Alert (Kraken)
  "Centralisation"         = "#8B0000", # Deep Crimson (Kraken)
  "Functional Redundancy"  = "#99D9D9" # Ice Blue (Kraken)
  # these are extra colours I might want to use...
  #"Shadow Blue"            = "#68A2B9"  # Shadow Blue (Kraken)
)

# Converts the named vector into a dataframe for easier mapping in ggplot
pal_df <- tibble(
  label   = names(kraken_palette), # printed label
  colour = kraken_palette, # colour code
  value   = factor(1:length(kraken_palette)) # module number
)

module_names <- tibble(
  Cluster = 1:7,
  Module_Name = c(
    "Macro Complexity",
    "Trophic Integration",
    "Energy Transport",
    "Trophic Asymmetry",
    "Control Heterogeneity",
    "Centralisation",
    "Functional Redundancy"
  )
)

# continuous ramp
seattle_anchors <- c(
  "Foam"      = "#E9E3D3", 
  "Ice"       = "#99D9D9", 
  "Shadow"    = "#68A2B9", 
  "Boundless" = "#355464", 
  "DeepSea"   = "#001628"
)

# 2. Create the Generator Function
# This creates a function that can interpolate any number of colors
seattle_abyssal_gen <- colorRampPalette(seattle_anchors)

# diverging colour
pal_diverge <- tibble(low = "#2C3E50",
                      mid = "white",
                      high = "#C1785A")

figure_theme <- function() {
  theme_bw() %+replace% 
    theme(
      # Text & Axes
      text = element_text(family = "space", color = "#001628"),
      plot.title = element_text(face = "bold", size = 14, margin = ggplot2::margin(b = 10)),
      axis.text = element_text(size = 10, color = "#1A1A1A"),
      axis.title = element_text(face = "bold"),
      
      # Grid & Background
      panel.grid.major = element_line(color = "#A5ACAF", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      
      # Border & Strip Alignment
      panel.border = element_rect(color = "#1A1A1A", fill = NA, linewidth = 1),
      strip.background = element_rect(color = "#1A1A1A", fill = "white", linewidth = 1),
      
      # Important: ensure there is no gap between strip and panel
      panel.spacing = unit(0, "pt"),
      
      # Legend
      legend.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}
