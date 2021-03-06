% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{get_map2}
\alias{get_map2}
\title{Alternative \code{get_map} function}
\usage{
get_map2(
  location = c(lon = -95.3632715, lat = 29.7632836),
  zoom = "auto",
  scale = "auto",
  maptype = c("terrain", "terrain-background", "satellite", "roadmap", "hybrid",
    "toner", "watercolor", "terrain-labels", "terrain-lines", "toner-2010", "toner-2011",
    "toner-background", "toner-hybrid", "toner-labels", "toner-lines", "toner-lite"),
  source = c("google", "osm", "stamen"),
  force = ifelse(source == "google", TRUE, FALSE),
  messaging = FALSE,
  urlonly = FALSE,
  filename = NULL,
  crop = TRUE,
  color = c("color", "bw"),
  language = "en-EN",
  ...
)
}
\arguments{
\item{location}{an address, longitude/latitude pair (in that order), or
left/bottom/right/top bounding box}

\item{zoom}{map zoom, an integer from 3 (continent) to 21 (building), default
value 10 (city).  openstreetmaps limits a zoom of 18, and the limit on
stamen maps depends on the maptype.  "auto" automatically determines the
zoom for bounding box specifications, and is defaulted to 10 with
center/zoom specifications.  maps of the whole world currently not
supported.}

\item{scale}{scale argument of get_googlemap() or get_openstreetmap()}

\item{maptype}{character string providing map theme. options available are
"terrain", "terrain-background", "satellite", "roadmap", and "hybrid"
(google maps), "terrain", "watercolor", and "toner" (stamen maps)}

\item{source}{Google Maps ("google"), OpenStreetMap ("osm"), Stamen Maps
("stamen")}

\item{force}{force new map (don't use archived version)}

\item{messaging}{turn messaging on/off}

\item{urlonly}{return url only}

\item{filename}{destination file for download (file extension added according
to format). Default \code{NULL} means a random tempfile().}

\item{crop}{(stamen and cloudmade maps) crop tiles to bounding box}

\item{color}{color ("color") or black-and-white ("bw")}

\item{language}{language for google maps}

\item{...}{...}
}
\value{
a ggmap object (a classed raster object with a bounding box
attribute)
}
\description{
An alternative to \code{ggmap}'s \code{get_map} function
}
\details{
An error in the CRAN available version of get_map means it will only plot Google Maps objects rather than Stamen or OSM maps. When ggmap is
updated this function will be removed. See \code{help(get_map)} for more details.
}
