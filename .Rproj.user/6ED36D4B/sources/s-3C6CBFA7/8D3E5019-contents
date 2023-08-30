

library(Rraven)

library(warbleR)

warbleR_options(wav.path = "/home/m/Dropbox/Projects/Phyllostomid call evolution/grabaciones", wl = 150, ovlp = 90, flim = c(30, 135))

sels <- imp_raven(path = .Options$warbleR$wav.path, sound.file.col = "Begin File", all.data = FALSE, waveform = FALSE)

sels <- sels[!duplicated(sels), ]

sa <- spec_an(sels)

spectrograms(sels, mar = mean(sa$duration)/2, pal = reverse.heat.colors)
trackfreqs(sels, length.out = 20, mar = mean(sa$duration)/2, bp = c(30, 135), ff.method = "tuneR")

st <- sel_tailor(sels, mar = mean(sa$duration)/2, frange = TRUE, auto.next = TRUE)
