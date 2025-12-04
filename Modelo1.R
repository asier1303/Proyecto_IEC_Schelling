# Cargar librerías necesarias
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("viridis")) install.packages("viridis")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("grid")) install.packages("grid")

library(ggplot2)
library(viridis)
library(gridExtra)
library(grid)

# ====================================================================================
# FUNCIÓN DE SIMULACIÓN
# ====================================================================================

simular_modelo1 <- function(N, p, citySize, years, alikePref, alfa){
  
  T <- citySize**2
  density <- round(N/T,3)   
  N1 <- round(N*p)
  N2 <- N-N1
  
  groups <- c(rep(1, N1), rep(2, N2))
  empty <- c(rep(0,T-N))
  set.seed(1234)
  city <- matrix(sample(c(groups, empty), citySize^2, replace = F), ncol = citySize)
  city0 <- city
  city0.des <- matrix(city0,T,1)
  
  getNeighbors <- function(coords) {
    n <- c()
    for (i in c(1:8)) {
      if (i == 1) { x <- coords[1] + 1; y <- coords[2] }
      if (i == 2) { x <- coords[1] + 1; y <- coords[2] + 1 }
      if (i == 3) { x <- coords[1]; y <- coords[2] + 1 }
      if (i == 4) { x <- coords[1] - 1; y <- coords[2] + 1 }
      if (i == 5) { x <- coords[1] - 1; y <- coords[2] }
      if (i == 6) { x <- coords[1] - 1; y <- coords[2] - 1 }
      if (i == 7) { x <- coords[1]; y <- coords[2] - 1 }
      if (i == 8) { x <- coords[1] + 1; y <- coords[2] - 1 }
      
      if (x < 1) x <- citySize
      if (x > citySize) x <- 1
      if (y < 1) y <- citySize
      if (y > citySize) y <- 1
      n <- rbind(n, c(x,y))
    }
    n
  }
  
  unHappyMonitor <- c()
  similarMonitor <- c()
  S <- c()
  I <- c()
  
  for (t in c(1:years)) {
    
    happy <- c()
    unhappy <- c() 
    similar <- c()
    similara <- c()
    similarb <- c()
    similar3 <- c()
    
    alikePref1 <- rep(alikePref, years)
    alikePref2 <- rep(alikePref*alfa, years)
    
    for (j in c(1:citySize)) {
      for (k in c(1:citySize)) {
        
        current <- c(j,k)
        value <- city[j,k] 
        
        if (value > 0) {
          
          likeNeighbors <- 0
          allNeighbors <- 0
          neighbors <- getNeighbors(current)
          
          for (i in c(1:nrow(neighbors))){
            x <- neighbors[i,1]
            y <- neighbors[i,2]
            
            if (city[x,y] > 0) {
              allNeighbors <- allNeighbors + 1
            }
            
            if (city[x,y] == value) {
              likeNeighbors <- likeNeighbors + 1
            }
          }
          
          if (is.nan(likeNeighbors / allNeighbors) == FALSE) {
            
            similar <- c(similar, likeNeighbors / allNeighbors)
            similara <- c(similara, likeNeighbors)
            similarb <- c(similarb, allNeighbors)
            
            if (value==1){
              alikePrefg <- alikePref1[t]
            }
            if (value==2){
              alikePrefg <- alikePref2[t]
            }
            
            similar1 <- likeNeighbors / N1
            similar2 <- (allNeighbors-likeNeighbors) / N2
            similar3 <- c(similar3, abs(similar1-similar2))
            
            if ((likeNeighbors / allNeighbors) < alikePrefg) {
              unhappy <- rbind(unhappy, current)
            } else {
              happy <- rbind(happy, current)
            }
            
          } else {
            happy <- rbind(happy, current)
          }
        }
      }
    }
    
    unHappyMonitor <- c(unHappyMonitor, length(unhappy)/(length(happy) + length(unhappy)))
    similarMonitor <- c(similarMonitor, mean(similar))
    S <- c(S, sum(similara)/sum(similarb))
    I <- c(I, 0.5*sum(similar3))
    
    if (length(unhappy) > 0) {
      randUnHappy <- sample(nrow(unhappy))
      
      for (i in randUnHappy) {
        mover <- unhappy[i,]
        mover_val <- city[mover[1],mover[2]]
        move_to <- c(sample(1:citySize,1), sample(1:citySize,1))
        move_to_val <- city[move_to[1], move_to[2]]
        
        while (move_to_val > 0) {
          move_to <- c(sample(1:citySize,1), sample(1:citySize,1))
          move_to_val <- city[move_to[1], move_to[2]]
        }
        
        city[mover[1], mover[2]] <- 0
        city[move_to[1], move_to[2]] <- mover_val
      }
    }
  }
  
  return(list(
    unhappy = unHappyMonitor, 
    similar = similarMonitor,
    S = S,
    I = I,
    city_final = city
  ))
}

# ====================================================================================
# FUNCIÓN PRINCIPAL DE VISUALIZACIÓN
# ====================================================================================

visualizar_schelling <- function(N = 2000, p = 0.5, citySize = 50, 
                                 years = 200, alikePref = 0.6, alfa = 1.5,
                                 guardar_graficos = FALSE) {
  
  # Iniciar temporizador total
  tiempo_inicio_total <- Sys.time()
  
  T <- citySize^2
  if (N > T) {
    stop(sprintf("Error: N (%d) no puede ser mayor que citySize^2 (%d)", N, T))
  }
  if (p < 0 || p > 1) {
    stop("Error: p debe estar entre 0 y 1")
  }
  if (alikePref < 0 || alikePref > 1) {
    stop("Error: alikePref debe estar entre 0 y 1")
  }
  
  cat("\n=== EJECUTANDO SIMULACIÓN ===\n")
  cat(sprintf("Parámetros: N=%d, p=%.2f, citySize=%d, years=%d, alikePref=%.2f, alfa=%.2f\n", 
              N, p, citySize, years, alikePref, alfa))
  
  # Medir tiempo de simulación
  tiempo_inicio_sim <- Sys.time()
  resultado <- simular_modelo1(N, p, citySize, years, alikePref, alfa)
  tiempo_fin_sim <- Sys.time()
  tiempo_simulacion <- as.numeric(difftime(tiempo_fin_sim, tiempo_inicio_sim, units = "secs"))
  
  cat(sprintf("✓ Simulación completada en %.2f segundos\n\n", tiempo_simulacion))
  
  # Medir tiempo de generación de gráficos
  tiempo_inicio_graficos <- Sys.time()
  
  tiempo <- 1:length(resultado$unhappy)
  diferencia <- resultado$unhappy - resultado$similar
  idx_cero <- which.min(abs(diferencia))
  
  # ====================================================================================
  # GRÁFICO 1: MAPA DE CALOR DE SEGREGACIÓN
  # ====================================================================================
  
  cat("📊 Generando Gráfico 1: Mapa de Calor...\n")
  
  city_df <- as.data.frame(as.table(resultado$city_final))
  names(city_df) <- c("X", "Y", "Grupo")
  city_df$Grupo <- factor(city_df$Grupo, 
                          levels = c(0, 1, 2),
                          labels = c("Vacío", "Grupo 1", "Grupo 2"))
  
  grafico1 <- ggplot(city_df, aes(x = Y, y = X, fill = Grupo)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_manual(values = c("Vacío" = "white", 
                                 "Grupo 1" = "red", 
                                 "Grupo 2" = "blue")) +
    coord_equal() +
    labs(title = sprintf("N = %d   p = %.1f   alpha = %.1f", N, p, alfa),
         x = "", y = "") +
    theme_minimal(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = "black", linewidth = 1),
      panel.background = element_rect(fill = "white", color = "black", linewidth = 2),
      panel.grid.minor = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "top",
      legend.title = element_blank()
    )
  
  print(grafico1)
  if (guardar_graficos) {
    ggsave("mapa_calor_segregacion.png", grafico1, width = 6, height = 6, dpi = 300)
  }
  
  # ====================================================================================
  # GRÁFICO 2: EVOLUCIÓN TEMPORAL - DOS CURVAS
  # ====================================================================================
  
  cat("📊 Generando Gráfico 2: Evolución Temporal...\n")
  
  df_temporal2 <- data.frame(
    Tiempo = rep(tiempo, 3),
    Valor = c(resultado$unhappy, resultado$similar, resultado$S),
    Variable = factor(rep(c("Inf", "S1", "S2"), each = length(tiempo)),
                      levels = c("Inf", "S1", "S2"))
  )
  
  cruces <- which(diff(sign(diferencia)) != 0)
  puntos_cruce <- data.frame()
  
  if (length(cruces) > 0) {
    for (idx in cruces) {
      t1 <- tiempo[idx]
      t2 <- tiempo[idx + 1]
      y1 <- diferencia[idx]
      y2 <- diferencia[idx + 1]
      
      t_cruce <- t1 - y1 * (t2 - t1) / (y2 - y1)
      valor_cruce <- resultado$unhappy[idx] + 
        (resultado$unhappy[idx + 1] - resultado$unhappy[idx]) * 
        (t_cruce - t1) / (t2 - t1)
      
      puntos_cruce <- rbind(puntos_cruce, 
                            data.frame(Tiempo = t_cruce, Valor = valor_cruce))
    }
  }
  
  grafico2 <- ggplot(df_temporal2, aes(x = Tiempo, y = Valor)) +
    geom_line(aes(color = Variable, linetype = Variable), linewidth = 1.2) +
    scale_color_manual(
      name = "",
      values = c("Inf" = "#E74C3C", "S1" = "#27AE60", "S2" = "#3498DB"),
      labels = c("Inf" = "Infelices", "S1" = "Similares (s1)", "S2" = "Similares (s2)")
    ) +
    scale_linetype_manual(
      name = "",
      values = c("Inf" = "solid", "S1" = "solid", "S2" = "dashed"),
      labels = c("Inf" = "Infelices", "S1" = "Similares (s1)", "S2" = "Similares (s2)")
    ) +
    labs(
      title = "Evolución Temporal: Infelices vs Similares",
      subtitle = sprintf("N=%d, p=%.1f, alikePref=%.1f, α=%.1f", N, p, alikePref, alfa),
      x = "Tiempo (años)",
      y = "Proporción"
    ) +
    ylim(0, 1) +
    theme_minimal(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "white"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40")
    )
  
  if (nrow(puntos_cruce) > 0) {
    grafico2 <- grafico2 +
      geom_point(data = puntos_cruce, aes(x = Tiempo, y = Valor),
                 color = "#9B59B6", size = 3, shape = 21, fill = "white", stroke = 2) +
      geom_vline(data = puntos_cruce, aes(xintercept = Tiempo),
                 linetype = "dotted", color = "#9B59B6", alpha = 0.5)
  }
  
  punto_min <- data.frame(
    Tiempo = tiempo[idx_cero],
    Valor_inf = resultado$unhappy[idx_cero],
    Valor_sim = resultado$similar[idx_cero]
  )
  
  grafico2 <- grafico2 +
    geom_point(data = punto_min, aes(x = Tiempo, y = Valor_inf),
               color = "#E74C3C", size = 3, shape = 19) +
    geom_point(data = punto_min, aes(x = Tiempo, y = Valor_sim),
               color = "#27AE60", size = 3, shape = 19) +
    geom_vline(xintercept = punto_min$Tiempo, 
               linetype = "dashed", color = "gray50", alpha = 0.7)
  
  print(grafico2)
  if (guardar_graficos) {
    ggsave("evolucion_temporal.png", grafico2, width = 10, height = 6, dpi = 300)
  }
  
  # ====================================================================================
  # GRÁFICO 3: ÍNDICE DE SEGREGACIÓN I 
  # ====================================================================================
  
  cat("📊 Generando Gráfico 3: Índice de Segregación...\n")
  
  I_normalizado <- resultado$I / (2 * max(resultado$I))
  
  df_segregacion <- data.frame(
    Tiempo = tiempo,
    I = I_normalizado
  )
  
  grafico3 <- ggplot(df_segregacion, aes(x = Tiempo, y = I)) +
    geom_line(color = "black", linewidth = 1) +
    labs(title = "",
         x = "Tiempo",
         y = "Segregación") +
    ylim(0, 1) +
    theme_minimal(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "white"),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = max(tiempo) * 0.6, y = 0.9,
             label = sprintf("Max I = %.2f", max(resultado$I)),
             hjust = 0, size = 4)
  
  print(grafico3)
  if (guardar_graficos) {
    ggsave("indice_segregacion.png", grafico3, width = 8, height = 5, dpi = 300)
  }
  
  # ====================================================================================
  # GRÁFICO 4: FUNCIÓN OBJETIVO
  # ====================================================================================
  
  cat("📊 Generando Gráfico 4: Función Objetivo...\n")
  
  df_objetivo <- data.frame(
    Tiempo = tiempo,
    Diferencia = diferencia
  )
  
  punto_cero <- data.frame(
    Tiempo = tiempo[idx_cero],
    Diferencia = diferencia[idx_cero]
  )
  
  grafico4 <- ggplot(df_objetivo, aes(x = Tiempo, y = Diferencia)) +
    geom_line(color = "#9B59B6", linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_point(data = punto_cero, aes(x = Tiempo, y = Diferencia),
               color = "#27AE60", size = 4) +
    geom_text(data = punto_cero, 
              aes(x = Tiempo, y = Diferencia, 
                  label = sprintf("(%.0f, %.4f)", Tiempo, Diferencia)),
              vjust = -1.5, hjust = 0.5, size = 3.5) +
    labs(title = "f(t) = Infelices - Similares",
         x = "Tiempo",
         y = "f(t)") +
    theme_minimal(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "white"),
      panel.grid.minor = element_blank()
    )
  
  print(grafico4)
  if (guardar_graficos) {
    ggsave("funcion_objetivo.png", grafico4, width = 8, height = 5, dpi = 300)
  }
  
  # ====================================================================================
  # GRÁFICO 5: PANEL COMPLETO (4 GRÁFICOS JUNTOS)
  # ====================================================================================
  
  cat("📊 Generando Gráfico 5: Panel completo (4 gráficos)...\n")
  
  # A) Mapa de calor simplificado
  grafico_panel_a <- ggplot(city_df, aes(x = Y, y = X, fill = Grupo)) +
    geom_tile(color = "white", linewidth = 0.05) +
    scale_fill_manual(values = c("Vacío" = "white", 
                                 "Grupo 1" = "red", 
                                 "Grupo 2" = "blue")) +
    coord_equal() +
    labs(title = "A) Distribución Espacial", x = "", y = "") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      panel.grid = element_blank()
    )
  
  # B) Evolución temporal
  grafico_panel_b <- ggplot(df_temporal2, aes(x = Tiempo, y = Valor)) +
    geom_line(aes(color = Variable, linetype = Variable), linewidth = 1) +
    scale_color_manual(
      name = "",
      values = c("Inf" = "#E74C3C", "S1" = "#27AE60", "S2" = "#3498DB"),
      labels = c("Inf" = "Infelices", "S1" = "s1", "S2" = "s2")
    ) +
    scale_linetype_manual(
      name = "",
      values = c("Inf" = "solid", "S1" = "solid", "S2" = "dashed"),
      labels = c("Inf" = "Infelices", "S1" = "s1", "S2" = "s2")
    ) +
    labs(title = "B) Evolución Temporal", x = "Tiempo (años)", y = "Proporción") +
    ylim(0, 1) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      legend.text = element_text(size = 7)
    )
  
  if (nrow(puntos_cruce) > 0) {
    grafico_panel_b <- grafico_panel_b +
      geom_vline(data = puntos_cruce, aes(xintercept = Tiempo),
                 linetype = "dotted", color = "#9B59B6", alpha = 0.5)
  }
  
  grafico_panel_b <- grafico_panel_b +
    geom_vline(xintercept = punto_min$Tiempo, 
               linetype = "dashed", color = "gray50", alpha = 0.7)
  
  # C) Índice de segregación
  grafico_panel_c <- ggplot(df_segregacion, aes(x = Tiempo, y = I)) +
    geom_line(color = "black", linewidth = 1) +
    labs(title = "C) Índice de Segregación I", x = "Tiempo (años)", y = "Segregación") +
    ylim(0, 1) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10)) +
    annotate("text", x = max(tiempo) * 0.6, y = 0.9,
             label = sprintf("Max I = %.2f", max(resultado$I)),
             hjust = 0, size = 3)
  
  # D) Función objetivo
  grafico_panel_d <- ggplot(df_objetivo, aes(x = Tiempo, y = Diferencia)) +
    geom_line(color = "#9B59B6", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_point(data = punto_cero, aes(x = Tiempo, y = Diferencia),
               color = "#27AE60", size = 2) +
    geom_vline(xintercept = punto_min$Tiempo, 
               linetype = "dashed", color = "gray50", alpha = 0.7) +
    geom_text(data = punto_cero, 
              aes(x = Tiempo, y = Diferencia, 
                  label = sprintf("(%.0f, %.4f)", Tiempo, Diferencia)),
              vjust = -1, hjust = 0.5, size = 2.5) +
    labs(title = "D)F(t) = Infelices - Similares", x = "Tiempo (años)", y = "f(t)") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10))
  
  # Crear panel combinado
  grafico5 <- grid.arrange(
    grafico_panel_a, grafico_panel_b, 
    grafico_panel_c, grafico_panel_d,
    ncol = 2
  )
  
  if (guardar_graficos) {
    ggsave("panel_completo_4_graficos.png", grafico5, width = 12, height = 10, dpi = 300)
  }
  
  # Finalizar temporizador de gráficos
  tiempo_fin_graficos <- Sys.time()
  tiempo_graficos <- as.numeric(difftime(tiempo_fin_graficos, tiempo_inicio_graficos, units = "secs"))
  
  # Tiempo total
  tiempo_fin_total <- Sys.time()
  tiempo_total <- as.numeric(difftime(tiempo_fin_total, tiempo_inicio_total, units = "secs"))
  
  # ====================================================================================
  # ESTADÍSTICAS RESUMEN
  # ====================================================================================
  
  cat("\n=== RESULTADOS ===\n\n")
  cat(sprintf("Parámetros:\n"))
  cat(sprintf("  N = %d agentes (%.1f%% ocupación)\n", N, 100*N/(citySize^2)))
  cat(sprintf("  Grupo 1: %.0f%% (tolerancia = %.2f)\n", 100*p, alikePref))
  cat(sprintf("  Grupo 2: %.0f%% (tolerancia = %.2f)\n", 100*(1-p), alikePref*alfa))
  cat(sprintf("  Ciudad: %dx%d parcelas\n", citySize, citySize))
  cat(sprintf("  Años simulados: %d\n\n", years))
  
  cat(sprintf("Resultados:\n"))
  cat(sprintf("  Tiempo óptimo (f(t) ≈ 0):              t = %d años\n", 
              tiempo[idx_cero]))
  cat(sprintf("  Valor de f(t*):                         %.6f\n", 
              diferencia[idx_cero]))
  cat(sprintf("  Índice I máximo:                        %.2f\n", max(resultado$I)))
  cat(sprintf("  Índice I final:                         %.2f\n", resultado$I[length(resultado$I)]))
  
  if (nrow(puntos_cruce) > 0) {
    cat(sprintf("\n  Puntos de cruce detectados:             %d\n", nrow(puntos_cruce)))
    for (i in 1:nrow(puntos_cruce)) {
      cat(sprintf("    Cruce %d: t = %.2f años, valor = %.4f\n", 
                  i, puntos_cruce$Tiempo[i], puntos_cruce$Valor[i]))
    }
  }
  
  cat(sprintf("\n  Proporción final de infelices:          %.4f (%.1f%%)\n", 
              resultado$unhappy[length(resultado$unhappy)],
              100*resultado$unhappy[length(resultado$unhappy)]))
  cat(sprintf("  Proporción final de similares:          %.4f (%.1f%%)\n", 
              resultado$similar[length(resultado$similar)],
              100*resultado$similar[length(resultado$similar)]))
  cat(sprintf("  Diferencia final:                       %.6f\n\n", 
              diferencia[length(diferencia)]))
  
  cat("✓ Todos los gráficos generados\n\n")
  
  # ====================================================================================
  # TIEMPOS DE EJECUCIÓN
  # ====================================================================================
  
  cat("=== TIEMPOS DE EJECUCIÓN ===\n\n")
  cat(sprintf("  Simulación:        %.2f segundos (%.1f%%)\n", 
              tiempo_simulacion, 100*tiempo_simulacion/tiempo_total))
  cat(sprintf("  Gráficos:          %.2f segundos (%.1f%%)\n", 
              tiempo_graficos, 100*tiempo_graficos/tiempo_total))
  cat(sprintf("  TOTAL:             %.2f segundos\n\n", tiempo_total))
  
  if (tiempo_total < 60) {
    cat(sprintf("⏱️  Tiempo total: %.1f segundos\n\n", tiempo_total))
  } else if (tiempo_total < 3600) {
    cat(sprintf("⏱️  Tiempo total: %.1f minutos (%.0f segundos)\n\n", 
                tiempo_total/60, tiempo_total))
  } else {
    cat(sprintf("⏱️  Tiempo total: %.2f horas (%.0f minutos)\n\n", 
                tiempo_total/3600, tiempo_total/60))
  }
  
  invisible(list(
    resultado = resultado,
    diferencia = diferencia,
    tiempo_optimo = tiempo[idx_cero],
    puntos_cruce = puntos_cruce,
    I_max = max(resultado$I),
    graficos = list(g1 = grafico1, g2 = grafico2, 
                    g3 = grafico3, g4 = grafico4, g5 = grafico5),
    tiempos = list(
      simulacion = tiempo_simulacion,
      graficos = tiempo_graficos,
      total = tiempo_total
    )
  ))
}

# ====================================================================================
# EJEMPLOS DE USO
# ====================================================================================

cat("\n" , paste(rep("=", 70), collapse = ""), "\n")
cat("MODELO 1 - CON MEDICIÓN DE TIEMPOS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Uso básico:\n")
cat("  visualizar_schelling()\n\n")

cat("Ejemplos con parámetros del paper:\n")
cat("  # Escenario 1: N=500, p=0.5, alfa=0.5\n")
cat("  visualizar_schelling(N=500, p=0.5, citySize=25, alikePref=0.6, alfa=0.5)\n\n")
cat("  # Escenario 3: N=1500, p=0.25, alfa=0.5\n")
cat("  visualizar_schelling(N=1500, p=0.25, citySize=40, alikePref=0.6, alfa=0.5)\n\n")

cat("Para guardar los gráficos:\n")
cat("  visualizar_schelling(guardar_graficos = TRUE)\n\n")

cat(paste(rep("=", 70), collapse = ""), "\n\n")

# EJECUTAR CON PARÁMETROS POR DEFECTO
visualizar_schelling()