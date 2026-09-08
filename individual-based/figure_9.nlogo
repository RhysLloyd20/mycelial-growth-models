; Mycelial growth model based on Hopkins & Boswell (2012)
; Line segments (turtles) grow, branch, and fuse (anastomosis) while resource is translocated through the network
; Tracks:
;   si        = internal resource
;   sext      = external resource
;   biomass   = hyphal segments per unit area
;   inhibitor = field biasing tips away from biomass

extensions [ rnd ]

breed [ segments segment ]

; Connection links segments that exchange resource (parent to child or anastomosis)
undirected-link-breed [ connections connection ]

connections-own [
  flux
]

globals [
  dt
  step-mm
  cell-mm
  per-mm
  step-len
  cell-area
  e-hat
  branch-thres
  b
  branch-scale
  De DC
  c1 c2 c3
  Tp Ta
  limiter-passes
  si0 se0
  si-hat
  turn-step
  turn-scale
  bias-coef
  d-hat
  sigma
  branch-mean
  branch-sd
  resource-rad
  supplement-rad
  n-spokes
  fuse-rad
  biomass-samples
  sample-offsets
]

patches-own [
  sext
  inhibitor
  biomass
  f0 f1
]

segments-own [
  seg-length
  si
  extended?
  branched?
  fused?
  par-link
  dist-tip
]

; Resource patch at centre with 8 line segments (one every 45 degrees)
to setup
  clear-all
  random-seed 123
  setup-world

  ; Time and space discretisation
  set dt 0.01              ; time step (d)
  set step-mm 0.05         ; segment length, dx (mm)
  set cell-mm 0.1          ; grid spacing, dy (mm) (Published value is 0.01 which is too slow to diffuse here)
  set per-mm (1 / cell-mm)
  set step-len (step-mm * per-mm)
  set cell-area (cell-mm * cell-mm)

  ; Growth and branching
  set e-hat 1e-12         ; minimum resource for extension (mol um^-2)
  set branch-thres 1e-11  ; minimum resource for branching (mol um^-2)
  set b 2.5e12            ; branching coefficient (given without unit)

  ; Eq. 9 as published is P = b * si * dt
  ; Returns 0.25 at branching threshold and saturates at 1 for si >~4e-11 meaning every eligible segment branches every step
  ; That doesn't produce the published networks
  ; branch-scale has been added as a scaling factor to reduce density
  set branch-scale (step-mm * step-mm)

  ; Diffusion coefficients (mm^2 d^-1)
  ; Table 1 gives 1 and 20 without units (below are equivalents at dy = 0.1)
  ; rk-diffuse converts to grid units by dividing by dy^2 (matching thesis Eq. 2.8)
  set De 0.01  ; external resource diffusion
  set DC 0.2   ; inhibitor diffusion

  ; Uptake (Michaelis-Menten (Eq. 3))
  set c1 600    ; uptake coefficient
  set c2 4e-9   ; uptake half-saturation
  set c3 0.01   ; inhibitor production rate

  ; Translocation (Eq. 1)
  set Tp 1                     ; passive/diffusive component
  set Ta active-translocation  ; active component (slider; Table 1 value is 20)
  set limiter-passes 5         ; passes used to hold fluxes within capacity

  ; Resource concentrations (mol um^-2)
  set si0 4e-9          ; initial internal resource
  set se0 3e-8          ; initial external resource
  set si-hat 6e-8       ; segment internal resource capacity

  ; Turning (velocity-jump (Eqs. 6-8))
  set d-hat sensitivity    ; sensitivity to inhibitor gradient (slider (Table 1 value is 20))
  set sigma 1.5            ; turning rate sd
  set turn-step 15         ; angular step size, dtheta = pi / 12 rad
  set turn-scale ((sigma * sigma * dt) / (2 * (pi / 12) * (pi / 12)))  ; total turn probability (must stay <1)
  set bias-coef (2 * d-hat / (sigma * sigma))

  ; Branch angle distribution (phi ~ N(0.44 pi, (pi / 18) ^ 2))
  set branch-mean (0.44 * 180)
  set branch-sd 10

  ; Initial geometry
  set resource-rad (0.5 * per-mm)
  set supplement-rad (0.3 * per-mm)
  set n-spokes 8
  set fuse-rad (2 * step-len + 1)   ; furthest crossing segment's end can be
  set biomass-samples 10

  ; Fractions along segment biomass gets sampled
  set sample-offsets n-values biomass-samples [ i -> (i + 0.5) / biomass-samples ]

  ; Inoculum resource
  ask patches with [ distancexy 0 0 <= resource-rad ] [ set sext se0 ]

  ; Central hub is a junction not a hypha
  let centre nobody
  create-segments 1 [
    new-segment 0 si0
    set seg-length 0
    set centre self
  ]

  ; Initial segments/spokes from centre
  ask centre [
    foreach (n-values n-spokes [ k -> k * 360 / n-spokes ]) [ h ->
      hatch-segments 1 [
        new-segment h si0
        forward step-len
        connect-to myself
      ]
    ]
  ]

  update-top
  update-biomass
  draw-network
  reset-ticks
end

to setup-world
  ask patches [ set pcolor white ]
end

; Order follows thesis appendix
to go
  if ticks >= max-ticks [ stop ]
  apply-supplement
  diffuse-external
  extend-tips
  branch-segments
  update-top
  update-biomass
  translocate
  uptake
  diffuse-inhibitor
  draw-network
  tick
end

; To add extra resource at set time
to apply-supplement
  if not add-supplement? [ stop ]
  if ticks != (round (supplement-day / dt)) [ stop ]

  ask patches with [ distancexy (supplement-x-mm * per-mm) (supplement-y-mm * per-mm) <= supplement-rad ] [
    set sext (sext + se0)
  ]
end

; Each tip with enough resource builds new segment
; Splits internal resource between old and new
; Paper only subtracts e-hat from parent but leaving child empty holds it below e-hat so it can't extend next step
; Halving the 0.5 cm/d growth rate the thesis calibrates to (Sec. 2.1.6)
to extend-tips
  ask segments with [ is-tip? and si > e-hat ] [
    if spawn-child new-heading [ set extended? true ]
  ]
end

; Sub-apical branching: segments with enough resource can spawn new segment at angle +/- phi (Eq. 9)
to branch-segments
  ask segments with [ can-branch? ] [
    if random-float 1 < branch-prob [
      let phi clamp-val (random-normal branch-mean branch-sd) 1 179
      let th (heading + (ifelse-value (random-float 1 < 0.5) [ 0 - phi ] [ phi ]))
      if spawn-child th [ set branched? true ]
    ]
  ]
end

to-report spawn-child [ th ]
  let remaining (si - e-hat)
  let child nobody

  hatch-segments 1 [
    new-segment th (remaining / 2)
    ifelse can-place?
      [ forward step-len  connect-to myself  set child self ]
      [ die ]
  ]

  if child = nobody [ report false ]

  set si (remaining / 2)
  ask child [ check-anastomosis myself ]
  report true
end

; Branching probability (Eq. 9)
to-report branch-prob
  report min (list 1 (b * branch-scale * si * dt))
end

; Initialise new segment where parent ends
to new-segment [ h sub ]
  set heading h
  set seg-length step-len
  set si sub
  set extended? false
  set branched? false
  set fused? false
  set par-link nobody
  hide-turtle
end

; Connect to parent
to connect-to [ par ]
  create-connection-with par
  set par-link connection-with par
end

; Anastomosis if new segment crosses existing one
; Obstructing tip is unaffected and continues (thesis Sec. 2.1.2, Fig. 2.3)
to check-anastomosis [ par ]
  let others (other segments in-radius fuse-rad) with [ self != par ]
  if not any? others [ stop ]

  let target min-one-of others [ crossing-of myself ]
  let t [ crossing-of myself ] of target
  if t > 1 [ stop ]

  jump (0 - (1 - t) * seg-length)
  set seg-length (t * seg-length)
  set fused? true
  create-connection-with target [ hide-link ]
end

; How far along segment s this segment crosses it (as fraction of the length of s) (2 if they miss)
; s is the extending segment so t is measured along it
to-report crossing-of [ s ]
  let ax [start-x] of s
  let ay [start-y] of s
  let rx ([xcor] of s - ax)
  let ry ([ycor] of s - ay)
  let qx (xcor - start-x)
  let qy (ycor - start-y)
  let denom (rx * qy - ry * qx)
  if abs denom < 1e-12 [ report 2 ]

  let px (start-x - ax)
  let py (start-y - ay)
  let t ((px * qy - py * qx) / denom)
  let u ((px * ry - py * rx) / denom)
  ifelse t > 1e-9 and t < (1 - 1e-9) and u > 1e-9 and u < (1 - 1e-9)
    [ report t ]
    [ report 2 ]
end

; Number each segment by how many links it is from the nearest tip (d(j) in Eq. 1)
; Translocation uses this to tell which way is towards a tip
; Segments that could still extend are the sinks (0)
to update-top
  ask segments [ set dist-tip false ]
  let frontier segments with [ is-tip? ]
  ask frontier [ set dist-tip 0 ]
  let d 0

  while [ any? frontier ] [
    set d (d + 1)
    let nxt (turtle-set [ connection-neighbors ] of frontier) with [ dist-tip = false ]
    ask nxt [ set dist-tip d ]
    set frontier nxt
  ]
end

; B in Eq. 5: shares each segment across every cell it crosses
to update-biomass
  let piece (1 / biomass-samples / cell-area)
  ask patches [ set biomass 0 ]
  ask segments with [ seg-length > 0 ] [
    foreach sample-offsets [ f ->
      let p patch-at (0 - f * seg-length * dx) (0 - f * seg-length * dy)
      if p != nobody [ ask p [ set biomass (biomass + piece) ] ]
    ]
  ]
end

; Internal translocation: passive (Tp) toward lower concentration and active (Ta) toward tips (Eqs. 1-2)
; Flux is held on the link and applied with opposite signs at each end so internal resource is conserved
to translocate
  ask connections [
    let sj [si] of end1
    let sk [si] of end2
    let dj [dist-tip] of end1
    let dk [dist-tip] of end2

    let diffusive (Tp * (sj - sk))
    let active 0
    if dk < dj [ set active (Ta * (dj - dk) * sj) ]
    if dj < dk [ set active (0 - (Ta * (dk - dj) * sk)) ]

    set flux (dt * (diffusive + active))
  ]

  repeat limiter-passes [ limit-fluxes ]

  ; Clamping here means limiter didn't converge and resource isn't conserved
  ask segments [
    set si clamp-val (si + net-flux) 0 si-hat
  ]
end

; Scale back inflows to segments that would overfill and outflows from ones that would empty
to limit-fluxes
  ask segments [
    let net net-flux
    ifelse net > 0 [
      let headroom (si-hat - si)
      if net > headroom [
        let factor (headroom / net)
        ask my-connections with [ inflow myself ] [ set flux (flux * factor) ]
      ]
    ] [
      if (0 - net) > si [
        let factor (si / (0 - net))
        ask my-connections with [ not inflow myself ] [ set flux (flux * factor) ]
      ]
    ]
  ]
end

to-report net-flux
  report sum [ ifelse-value (end2 = myself) [ flux ] [ 0 - flux ] ] of my-connections
end

to-report inflow [ seg ]
  report ifelse-value (end2 = seg) [ flux > 0 ] [ flux < 0 ]
end

; Michaelis-Menten uptake from cell at start of the segment (Eq. 3)
; Eq. 3 says the start of the segment thesis Sec. 2.1.4 says the end (Following paper here)
; Oldest first (thesis Sec. 2.1.4)
; Eq. 3 as published is take = c1 * m * sext * dt
; c1 * m * dt reaches 3-6 at the calibrated c1, so that form removes several times the substrate present
; Exponential below is the solution of Eq. 4 reaction term over the step and cannot overdraw
to uptake
  foreach (sort (segments with [ seg-length > 0 ])) [ s ->
    ask s [
      let cell start-patch
      if cell != nobody [
        let avail [sext] of cell
        if avail > 0 [
          let m (si / (c2 + si))
          let take (min (list (avail * (1 - exp (0 - c1 * m * dt)))
                              (si-hat - si)))
          if take > 0 [
            set si (si + take)
            ask cell [ set sext (sext - take) ]
          ]
        ]
      ]
    ]
  ]
end

; External resource diffusion (Eq. 4)
to diffuse-external
  ask patches [ set f0 sext ]
  rk-diffuse De
  ask patches [ set sext f0 ]
end

; Inhibitor production (proportional to local biomass) and diffusion (Eq. 5)
to diffuse-inhibitor
  ask patches [ set inhibitor (inhibitor + dt * c3 * biomass) ]
  ask patches [ set f0 inhibitor ]
  rk-diffuse DC
  ask patches [ set inhibitor f0 ]
end

; Second-order Runge-Kutta diffusion step using discrete Laplacian
; Coefficient divided by dy^2 (thesis Eq. 2.8)
; Sub-divided to stay inside stability limit
to rk-diffuse [ Dcoef ]
  if Dcoef <= 0 [ stop ]
  let k (Dcoef / cell-area)
  let n ceiling (8 * k * dt)
  let h (dt / n)
  repeat n [
    ask patches [ set f1 (f0 + 0.5 * h * k * lap-f0) ]
    ask patches [ set f0 max (list 0 (f0 + h * k * lap-f1)) ]
  ]
end

; Laplacian in grid cells
; Zero-flux boundaries applied by copying missing neighbour (thesis Eqs. 2.9-2.10)
; Interior patches take fast path so only the edges need copying
to-report lap-f0
  if count neighbors4 = 4 [ report (sum [f0] of neighbors4) - 4 * f0 ]
  report (mirror-f0 1 0) + (mirror-f0 -1 0) + (mirror-f0 0 1) + (mirror-f0 0 -1) - 4 * f0
end

to-report lap-f1
  if count neighbors4 = 4 [ report (sum [f1] of neighbors4) - 4 * f1 ]
  report (mirror-f1 1 0) + (mirror-f1 -1 0) + (mirror-f1 0 1) + (mirror-f1 0 -1) - 4 * f1
end

; Value of neighbour at this offset or of the opposite one where neighbour is off grid
to-report mirror-f0 [ ox oy ]
  let p patch-at ox oy
  if p != nobody [ report [f0] of p ]
  let q patch-at (0 - ox) (0 - oy)
  report ifelse-value (q = nobody) [ f0 ] [ [f0] of q ]
end

to-report mirror-f1 [ ox oy ]
  let p patch-at ox oy
  if p != nobody [ report [f1] of p ]
  let q patch-at (0 - ox) (0 - oy)
  report ifelse-value (q = nobody) [ f1 ] [ [f1] of q ]
end

; Colour and thickness of segments by internal resource concentration
; Log-scaled from e-hat to si-hat, green-orange-red (so frames are comparable)
to draw-network
  ask segments with [ par-link != nobody ] [
    let fr clamp-val ((ln (si + 1e-30) - ln e-hat) / (ln si-hat - ln e-hat)) 0 1
    let shade 0

    ; Green (low) to orange (med) to red (high)
    ifelse fr < 0.5 [
      ; Green to orange: red rises, green stays high
      let g (fr / 0.5)
      set shade (list (g * 255) (180 + (1 - g) * 75) 0)
    ] [
      ; Orange to red: green drops
      let g ((fr - 0.5) / 0.5)
      set shade (list 255 (180 * (1 - g)) 0)
    ]

    ask par-link [
      set color shade
      set thickness (0.05 + 0.15 * fr)
    ]
  ]
end

; Tip = segment that can still extend (not already extended or fused)
; Hub excluded
to-report is-tip?
  report (seg-length > 0) and (not extended?) and (not fused?)
end

; Segment can only branch once (thesis Sec. 2.1.2) with enough resource and only after extending
to-report can-branch?
  report (seg-length > 0) and (si >= branch-thres) and extended? and (not branched?)
end

; Tips stop at world edge
to-report can-place?
  report patch-ahead step-len != nobody
end

; Start of segment is 1 length back
to-report start-x
  report xcor - seg-length * dx
end

to-report start-y
  report ycor - seg-length * dy
end

to-report start-patch
  report patch-at (0 - seg-length * dx) (0 - seg-length * dy)
end

; Velocity-jump turning probabilities (Eqs. 7-8)
; Based on local inhibitor gradient biasing tips away from high inhibitor (dense growth)
to-report new-heading
  let gx (((inhib-at 1 0) - (inhib-at -1 0)) / 2)
  let gy (((inhib-at 0 1) - (inhib-at 0 -1)) / 2)
  let p-right (turn-scale / 2)
  let p-left  (turn-scale / 2)

  if (sqrt (gx * gx + gy * gy)) > 1e-12 [
    let away atan (0 - gx) (0 - gy)
    let tau-right tau-of (heading + turn-step / 2) away
    let tau-left  tau-of (heading - turn-step / 2) away
    let denom (tau-right + tau-left)
    set p-right (turn-scale * tau-right / denom)
    set p-left  (turn-scale * tau-left / denom)
  ]

  ; Turn clockwise, anti-clockwise or keep current direction
  let routes (list (list (heading + turn-step) p-right)
                   (list (heading - turn-step) p-left)
                   (list heading (1 - turn-scale)))

  report first rnd:weighted-one-of-list routes [ r -> last r ]
end

; Turning kernel (Eq. 8): exp(2 d-hat / sigma^2 * cos(phi - theta_p))
to-report tau-of [ h away ]
  report exp (bias-coef * cos (h - away))
end

to-report inhib-at [ ox oy ]
  let p patch-at ox oy
  report ifelse-value (p = nobody) [ inhibitor ] [ [inhibitor] of p ]
end

to-report clamp-val [ v lo hi ]
  report max (list lo (min (list hi v)))
end
@#$#@#$#@
GRAPHICS-WINDOW
312
23
1399
1111
-1
-1
10.6832
1
10
1
1
1
0
1
1
1
-50
50
-50
50
0
0
1
ticks
30.0

BUTTON
30
30
97
63
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
124
30
187
63
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
18
110
246
143
max-ticks
max-ticks
0
600
150.0
10
1
ticks
HORIZONTAL

SLIDER
18
163
245
196
active-translocation
active-translocation
0
50
20.0
1
1
Ta
HORIZONTAL

SLIDER
18
217
245
250
sensitivity
sensitivity
0
40
20.0
1
1
d-hat
HORIZONTAL

SWITCH
15
270
246
303
add-supplement?
add-supplement?
1
1
-1000

SLIDER
13
327
246
360
supplement-day
supplement-day
0
4
0.5
0.5
1
day
HORIZONTAL

SLIDER
12
384
249
417
supplement-x-mm
supplement-x-mm
-3
3
1.0
0.1
1
mm
HORIZONTAL

SLIDER
11
442
249
475
supplement-y-mm
supplement-y-mm
-3
3
2.0
0.1
1
mm
HORIZONTAL

MONITOR
22
513
79
558
day
ticks * dt
2
1
11

MONITOR
114
516
180
561
segments
count segments
0
1
11

MONITOR
23
576
80
621
tips
count segments with [is-tip?]
0
1
11

MONITOR
112
575
192
620
max Si / cap
(max [si] of segments) / si-hat
3
1
11

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.4.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
