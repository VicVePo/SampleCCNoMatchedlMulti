#' Cálculo de Tamaño Muestral para Estudios de Casos y Controles
#' 
#' @param alfa Nivel de significación
#' @param potencia Potencia deseada (opcional si se proporciona n_casos)
#' @param n_casos Número de casos (opcional si se proporciona potencia)
#' @param p0 Proporción de exposición en controles (requerido si metodo = "rho")
#' @param OR Odds Ratio a detectar
#' @param IC_superior Límite superior del intervalo de confianza del OR
#' @param IC_inferior Límite inferior del intervalo de confianza del OR (opcional si se proporciona EE)
#' @param EE Error estándar del coeficiente (opcional si se proporcionan los IC)
#' @param n_previo Tamaño de muestra del estudio previo (requerido para método "ee")
#' @param r Razón de controles a casos, por defecto es 1
#' @param metodo Método de cálculo: "rho" o "ee" (error estándar)
#' @param valores_rho Vector de valores rho para calcular tamaños muestrales (si metodo = "rho")
#' @return Un data frame con resultados según el método elegido
#' @export
SampleCCNoMatchedlMulti <- function(alfa, potencia = NULL, n_casos = NULL,
                                    p0 = NULL, OR,
                                    IC_superior = NULL, IC_inferior = NULL,
                                    EE = NULL, n_previo = NULL,
                                    r = 1, 
                                    metodo = "rho",
                                    valores_rho = seq(0, 0.9, by = 0.1)) {
  
  # Verificaciones según el método
  if (metodo == "rho") {
    if (is.null(p0)) {
      stop("Para método 'rho' se requiere p0")
    }
  } else if (metodo == "ee") {
    if (is.null(EE) && (is.null(IC_superior) || is.null(IC_inferior))) {
      stop("Para método 'ee' se requiere EE o ambos límites del IC")
    }
    if (is.null(n_previo)) {
      stop("Para método 'ee' se requiere n_previo")
    }
  } else {
    stop("El método debe ser 'rho' o 'ee'")
  }
  
  # Si se proporcionaron los IC, calcular el EE
  if (!is.null(IC_superior) && !is.null(IC_inferior)) {
    EE <- (log(IC_superior) - log(OR)) / qnorm(0.975)
    cat("Error Estándar calculado:", EE, "\n")
  }
  
  if (metodo == "rho") {
    p1 <- (OR * p0) / (1 - p0 + (OR * p0))
    
    z_alfa <- qnorm(1 - alfa/2)
    if (!is.null(potencia)) {
      z_beta <- qnorm(potencia)
      
      calcular_n <- function(rho) {
        numerador <- (z_alfa * sqrt((1 + 1/r) * p0 * (1 - p0)) + 
                       z_beta * sqrt(p1 * (1 - p1) + p0 * (1 - p0) / r))^2
        denominador <- (p1 - p0)^2 * (1 - rho^2)
        
        n_casos <- ceiling(numerador / denominador)
        n_controles <- ceiling(n_casos * r)
        
        return(c(n_casos, n_controles))
      }
      
      resultados <- data.frame(
        rho = valores_rho,
        n_casos = sapply(valores_rho, function(rho) calcular_n(rho)[1]),
        n_controles = sapply(valores_rho, function(rho) calcular_n(rho)[2])
      )
      
      resultados$tamaño_total <- resultados$n_casos + resultados$n_controles
      
    } else {
      # Cálculo de potencia con método rho
      calcular_potencia <- function(rho) {
        numerador <- (p1 - p0)^2 * (1 - rho^2)
        denominador <- sqrt((1 + 1/r) * p0 * (1 - p0))
        z_beta <- sqrt(n_casos * numerador) / denominador - z_alfa
        return(pnorm(z_beta))
      }
      
      resultados <- data.frame(
        rho = valores_rho,
        potencia = sapply(valores_rho, calcular_potencia)
      )
    }
  } else {
    # Método basado en error estándar
    z_alfa <- qnorm(1 - alfa/2)
    if (!is.null(potencia)) {
      z_gamma <- qnorm(potencia)
      n <- ceiling((z_alfa + z_gamma)^2 * n_previo * EE^2 / (log(OR))^2)
      n_casos <- ceiling(n / (1 + r))
      n_controles <- ceiling(n_casos * r)
      resultados <- data.frame(
        n_casos = n_casos,
        n_controles = n_controles,
        tamaño_total = n_casos + n_controles
      )
    } else {
      z_gamma <- sqrt(n_casos * (log(OR))^2 / (n_previo * EE^2)) - z_alfa
      potencia <- pnorm(z_gamma)
      resultados <- data.frame(
        potencia = potencia,
        n_casos = n_casos,
        n_controles = n_casos * r,
        tamaño_total = n_casos * (1 + r)
      )
    }
  }
  
  return(resultados)
}

#' Logística de Estudio de Casos y Controles
#'
#' @param n_casos Número de casos requeridos
#' @param n_controles Número de controles requeridos
#' @param tasa_identificacion_casos Tasa de identificación de casos por día
#' @param tasa_seleccion_controles Tasa de selección de controles por día
#' @param tasa_rechazo_casos Tasa de rechazo esperada para casos
#' @param tasa_rechazo_controles Tasa de rechazo esperada para controles
#' @param dias_laborables_mes Número de días laborables por mes
#' @return Una lista con los cálculos logísticos del estudio
#' @export
logistica_estudio_caso_control <- function(n_casos, 
                                           n_controles,
                                           tasa_identificacion_casos,
                                           tasa_seleccion_controles,
                                           tasa_rechazo_casos,
                                           tasa_rechazo_controles,
                                           dias_laborables_mes) {
  
  casos_a_contactar <- n_casos / (1 - tasa_rechazo_casos)
  controles_a_contactar <- n_controles / (1 - tasa_rechazo_controles)
  
  dias_para_casos <- ceiling(casos_a_contactar / tasa_identificacion_casos)
  dias_para_controles <- ceiling(controles_a_contactar / tasa_seleccion_controles)
  
  total_dias <- max(dias_para_casos, dias_para_controles)
  duracion_meses <- total_dias / dias_laborables_mes
  
  resultados <- list(
    casos_requeridos = n_casos,
    controles_requeridos = n_controles,
    casos_a_contactar = ceiling(casos_a_contactar),
    controles_a_contactar = ceiling(controles_a_contactar),
    dias_estimados = total_dias,
    meses_estimados = round(duracion_meses, 2)
  )
  
  # Imprimir resumen
  cat("
Resumen de logística del estudio:
")
  cat("----------------------------------------
")
  cat("Casos requeridos:", n_casos, "
")
  cat("Controles requeridos:", n_controles, "
")
  cat("Casos a contactar:", ceiling(casos_a_contactar),
      "(", tasa_rechazo_casos*100, "% de rechazo)
")
  cat("Controles a contactar:", ceiling(controles_a_contactar),
      "(", tasa_rechazo_controles*100, "% de rechazo)
")
  cat("Días necesarios:", total_dias, "
")
  cat("Meses necesarios:", round(duracion_meses, 2),
      "(", dias_laborables_mes, "días laborables por mes)
")
  cat("----------------------------------------
")
  
  return(invisible(resultados))
}
