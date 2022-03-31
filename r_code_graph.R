library(tidyverse)

PATH <- "/Users/danil/BSE/Term2/Networks/Final/networks_project/"
setwd(PATH)


model_1 <- read.csv("model1_data.csv")
model_2 <- read.csv("model2_data.csv")
model_3 <- read.csv("model3_data.csv")
model_4 <- read.csv("model4_data.csv")



avg_infection_length <- function(model) {
  infection_length <- model %>%
    group_by(p_infection,infection_period,recovery_period,percolation, simulation_number) %>%
    summarise(length_of_pandemic = max(time_period)) %>%
    ungroup() %>%
    group_by(p_infection,infection_period,recovery_period,percolation) %>%
    summarise(length_of_pandemic = mean(length_of_pandemic))
}

infection_length <- function(model){
  len <- model %>%
    group_by(p_infection, infection_period, recovery_period, percolation, simulation_number) %>%
    summarise(length_of_pandemic = max(time_period))
  
  len
}

plot_len_hist <- function(len, color = len$p_infection, n = 1, labs = "Probability of Infection"){
  
  true_names <- as.character(unique(color))
  values = c("blue", "red", "green", "purple")
  
  ggplot(len, aes(x = length_of_pandemic, fill = as.factor(color))) +
    geom_histogram(bins = 60) +
    theme_classic() +
    ggtitle(paste0("Model " ,n)) +
    scale_color_manual(labels = true_names, values = values[1:length(true_names)]) +
    xlab("Length of Disease") + 
    ylab("Number of simulations") + 
    labs( fill = labs) +
    theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10),
          axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 10),
          plot.title = element_text(size = 12))
  
  ggsave(
    filename = paste0(labs, n, ".png"),
    dpi = 300,
    plot = last_plot()
  )
  
}


len_1 <- infection_length(model_1)
len_3 <- infection_length(model_3)


infection_model_1 <- avg_infection_length(model_1)
infection_model_2 <- avg_infection_length(model_2)
infection_model_3 <- avg_infection_length(model_3)
infection_model_4 <- avg_infection_length(model_4)

plot_infection <- function(model, n = 1, sim_num = 0, j = 1){
  pivoted <- model %>%
    pivot_longer(c(susceptible_pct, infected_pct, recovered_pct), names_to = "Percent")
  
  models <- pivoted %>%
    group_by(parameter_combination)%>%
    group_split()
  
  ggplot(models[[n]][models[[n]]$simulation_number==sim_num,], aes(x = time_period, y = value, col = Percent)) +
    geom_line() +
    theme_classic() +
    ggtitle(paste0("Probability of infection = ", models[[n]]$p_infection[1],
                   " Infection period: ", models[[n]]$infection_period[1])) +
    scale_color_manual(labels = c("Infected", "Recovered", "Susceptible"), values = c("blue", "red", "green")) +
    xlab("Time") + 
    ylab("Proportion") + 
    labs(col = "Proportions",
         subtitle = paste0(" Recovery period: ",models[[n]]$recovery_period[1],
                           " Percolation: ", models[[n]]$percolation[1])) +
    theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10),
          axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 10),
          plot.title = element_text(size = 12))
  
  
  ggsave(
    filename = paste0("P_inf_",models[[n]]$p_infection[1] ,
                      "inf_per_",models[[n]]$infection_period[1],
                      "rec_per_", models[[n]]$recovery_period[1],
                      "perc_", models[[n]]$percolation[1],
                      "model", j,
                      ".png"),
    dpi = 300,
    plot = last_plot()
 )
    
}


model_danger <- function(model){
  model_1_dang <- model %>% 
    group_by(p_infection,infection_period,recovery_period,percolation, simulation_number) %>%
    summarise(min_suceptible = min(susceptible_pct),
              length_of_pandemic = max(time_period),
              dangerous = as.numeric(min(susceptible_pct) < 0.7)) 
  model_1_dang
}


model1_dang <- model_danger(model_1)
model2_dang <- model_danger(model_2)
model3_dang <- model_danger(model_3)
model4_dang <- model_danger(model_4)

ggsave(
  filename = paste0("Probability", 4, ".png"),
  dpi = 300,
  plot = last_plot()
)

plot_infection(model_1, n = 1, sim_num = 8, j = 1)

plot_infection(model_1, n = 39, sim_num = 17, j = 1)

plot_infection(model_1, n = 71, sim_num = 14, j = 1)

plot_infection(model_2, n = 1, sim_num = 10, j = 2)

plot_infection(model_2, n = 10, sim_num = 5, j = 2)

plot_infection(model_2, n = 11, sim_num = 18, j = 2)

plot_infection(model_2, n = 11, sim_num = 20, j = 2)

plot_infection(model_2, n = 12, sim_num = 20, j = 2)

plot_infection(model_2, n = 46, sim_num = 29, j = 2)

plot_infection(model_2, n = 47, sim_num = 10, j = 2)

plot_infection(model_2, n = 61, sim_num = 10, j = 2)

plot_infection(model_2, n = 62, sim_num = 10, j = 2)


plot_infection(model_3, n = 10, sim_num = 15, j = 3)

plot_infection(model_3, n = 20, sim_num = 15, j = 3)
plot_infection(model_3, n = 30, sim_num = 15, j = 3)

plot_infection(model_3, n = 49, sim_num = 15, j = 3)

plot_infection(model_4, n = 4, sim_num = 11, j = 4)

plot_infection(model_4, n = 22, sim_num = 11, j = 4)

plot_infection(model_4, n = 30, sim_num = 11, j = 4)

plot_infection(model_4, n = 34, sim_num = 11, j = 4)

plot_infection(model_4, n = 45, sim_num = 11, j = 4)

plot_infection(model_4, n = 58, sim_num = 11, j = 4)



