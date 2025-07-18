# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

heatindex_vec <- function(T, rh) {
    .Call(`_heatindex_heatindex_vec`, T, rh)
}

wetbulb_vec <- function(p, T, rh, psychrometric, icebulb, verbose, Lewis) {
    .Call(`_heatindex_wetbulb_vec`, p, T, rh, psychrometric, icebulb, verbose, Lewis)
}

rh_from_wetbulb_vec <- function(p, T, Tw, psychrometric, icebulb, verbose, Lewis) {
    .Call(`_heatindex_rh_from_wetbulb_vec`, p, T, Tw, psychrometric, icebulb, verbose, Lewis)
}

pvstarl <- function(T) {
    .Call(`_heatindex_pvstarl`, T)
}

pvstars <- function(T) {
    .Call(`_heatindex_pvstars`, T)
}

pvstar <- function(T) {
    .Call(`_heatindex_pvstar`, T)
}

qvstarl <- function(p, T) {
    .Call(`_heatindex_qvstarl`, p, T)
}

qvstars <- function(p, T) {
    .Call(`_heatindex_qvstars`, p, T)
}

qvstar <- function(p, T) {
    .Call(`_heatindex_qvstar`, p, T)
}

Tstarl <- function(pv) {
    .Call(`_heatindex_Tstarl`, pv)
}

Tstars <- function(pv) {
    .Call(`_heatindex_Tstars`, pv)
}

Tstar <- function(pv) {
    .Call(`_heatindex_Tstar`, pv)
}

