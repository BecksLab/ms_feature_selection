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
  "Primary"   = "#001628", # Deep Sea Blue
  "Secondary" = "#99D9D9", # Ice Blue
  "Tertiary"  = "#68A2B9", # Shadow Blue
  "Highlight" = "#E9072B"  # Red Alert
)

# Converts the named vector into a dataframe for easier mapping in ggplot
pal_df <- data.frame(l = names(kraken_palette), c = kraken_palette)

kraken_7 <- c(
  "Kraken Deep"    = "#001628", # Deep Sea Blue (Kraken)
  "Action Green"   = "#69BE28", # Action Green (Seahawks)
  "Storm Gold"     = "#FBD800", # Storm Gold (Storm)
  "Torrent Teal"   = "#005A5B", # Torrent Teal (Torrent)
  "Kraken Red"     = "#E9072B", # Red Alert (Kraken)
  "Deep Crimson"   = "#8B0000", # Deep Crimson (Kraken)
  "Kraken Ice"     = "#99D9D9", # Ice Blue (Kraken)
  "Shadow Blue"    = "#68A2B9"  # Shadow Blue (Kraken)
  
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
