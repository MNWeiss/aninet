#' Contacts between southern resident killer whales
#'
#' A dataset containing the number of contacts observed between members of J pod
#' Data were collected from an unmanned aerial vehicle over 11 hours of observation in 2019
#'
#' @format A symmetric square matrix, with rows and columns representing individuals, and the entried indicating the contact rates between them.
"srkw_contact"

#' Estimated maternal kinship between southern resident killer whales
#'
#' Pairwise relatedness coefficients between individual southern resident killer whales in J pod, based on maternal pedigrees.
#'
#' @format A symmetric square matrix containing pairwise relatedness estimates
"srkw_kinship"

#' Individual attributes of southern resident killer whales
#'
#' Contains individual information about the members of J pod in 2019
#'
#' @format A data.frame containing individual ID, age, sex, age-sex class, and matriline identity.
"srkw_attributes"

#' Sightings of southern resident killer whales over 10 days
#'
#' Sightings history of J pod over a 10 day observation period
#'
#' @format A day-by-individual matrix, with 1 indicating an individual was seen on a given day, and a 0 indicating they were not.
"srkw_sightings"

#' Synchronous surfacings between southern resident killer whales
#'
#' Number of synchronous surfacing events between members of J pod in the summer of 2019.
#'
#' @format A symmetric square matrix containing the observed number of synchronous surfacing events between each dyad.
"srkw_surfacing"

#' Dyadic sampling effort of southern resident killer whales
#'
#' Total observation time (in seconds) for each pair of whales in J pod during which interactions were observed.
#'
#' @format A symmetric square matrix containing the dyadic sampling efforts
"srkw_sampling"

#' Dyadic codetections of southern resident killer whales
#'
#' Portion of sampling time in which each pair of whales in J pod was observed together. A value of 1 would indicate whales were always seen together, while 0 would indicate they were never seen together.
#'
#' @format A symmetric square matrix containing the co-detection measures.
"srkw_codetection"

