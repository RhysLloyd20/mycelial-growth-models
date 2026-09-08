; Mycelial growth model based on Boswell et al. (2007)
; Cells track resource concentrations while links between cells are active/inactive hyphae
; Growing tips (turtles) have a biased random walk, branch and leave behind trails of new hyphae
;   si      = internal resource
;   sext    = external resource
;   ahyphae = active hyphae
;   ihyphae = inactive hyphae

extensions [ rnd ]

breed [ cells cell ]
breed [ tips tip ]

; Active hyphae (translocate resource)
; Inactive hyphae (awaiting degradation)
undirected-link-breed [ ahyphae ahypha ]
undirected-link-breed [ ihyphae ihypha ]

ahyphae-own [
  gdir   ; direction of tip when hypha formed
]

globals [
  delta-x v Dp b Di De c1 c2 c3 c4 c5 omega Da-rich Da-poor d-i
  si0 se0
  t-end max-ticks target-pmove
  win-x0 win-x1 win-y0 win-y1
  cen-x0 cen-x1 cen-y0 cen-y1
  sec-x0 sec-x1 sec-y0 sec-y1
  dirs sqrt3-2
  simtime dt
]

cells-own [
  si sext              ; internal and external resource
  blk                  ; resource block index (-1 = none)
  neighbours           ; agentset of the 6 surrounding cells
  delta-si delta-se    ; accumulators for resource changes
  tip-here?
]

tips-own [
  tcell      ; cell this tip sits on
]

to setup
  clear-all
  random-seed 123
  set-default-params
  set sqrt3-2 (sqrt 3 / 2)

  ; Six hexagonal directions (rows run horizontally so neighbours are 30 degrees off vertical)
  set dirs [30 90 150 210 270 330]

  ; Lattice bounds (inset a cell from world edge)
  set win-x0 1
  set win-x1 53
  set win-y0 1
  set win-y1 53

  ; Resource blocks placement
  set cen-x0 15
  set cen-x1 23
  set cen-y0 12
  set cen-y1 20
  set sec-x0 31
  set sec-x1 39
  set sec-y0 34
  set sec-y1 42

  set simtime 0

  setup-world
  create-cell-grid
  find-neighbours
  setup-blocks
  setup-inoculum
  reset-ticks
end

; Parameters re-tuned. Boswell et al. (2007) give no values to reproduce their Fig. 6 (Only the c1:c3 ratio is from their paper)
; Lengths in mm, time in days, resource in model units
; Coefficients are converted to grid units at the point of use by dividing by delta-x
to set-default-params
  set delta-x    0.1       ; cell diameter (mm)
  set v          2.8       ; tip growth rate
  set Dp         0.004     ; tip diffusion
  set b          0.35      ; branch rate
  set Di         0.04      ; internal resource diffusion
  ; Glucose in agar is 34.56 mm^2 d^-1 (Olsson 1995). That value would set dt and make runs ~9x longer
  set De         0.003456  ; external resource diffusion
  set c1         8.0       ; uptake rate (external to internal)
  set c2         0.01      ; growth cost per unit hypha
  set c3        80.0       ; uptake cost (external depletion)
  set c4         0.001     ; active translocation cost
  set c5         0.8       ; maintenance cost
  set omega      1.0E-6    ; minimum resource for active hyphae
  set Da-rich    0.4       ; active translocation (resource-rich)
  set Da-poor    4.0       ; active translocation (resource-poor)
  set d-i        0.04      ; inactive hypha degradation rate
  set si0        1.5       ; initial internal resource
  set se0       70.0       ; initial external resource in blocks

  ; Run control
  set t-end        1.5     ; stop time (matches Fig. 6d)
  set max-ticks    300000  ; cap in case t-end is never reached
  set target-pmove 0.4     ; max tip movement probability per step
end

to setup-world
  resize-world 0 54 0 54
  set-patch-size 8
  ask patches [ set pcolor white ]
end

; Build hexagonal lattice: alternate rows are offset by half a cell width giving 60 degree branch angles
to create-cell-grid
  foreach (range (floor ((win-y1 - win-y0) / sqrt3-2) + 1)) [ i ->
    let row-y win-y0 + i * sqrt3-2
    let row-x win-x0 + 0.5 * (i mod 2)

    foreach (range (floor (win-x1 - row-x) + 1)) [ j ->
      create-cells 1 [
        setxy (row-x + j) row-y
        set blk -1
        set tip-here? false
        set shape "square"
        set size 0
      ]
    ]
  ]
end

; Each cell sits 1 lattice spacing from its six neighbours
to find-neighbours
  ask cells [ set neighbours other cells in-radius 1.01 ]
end

; Place square resource blocks on lattice
to setup-blocks
  let blocks (list (list cen-x0 cen-x1 cen-y0 cen-y1)
                   (list sec-x0 sec-x1 sec-y0 sec-y1))

  foreach (range length blocks) [ i ->
    let bd item i blocks
    ask cells with [ xcor >= item 0 bd and xcor <= item 1 bd and
                     ycor >= item 2 bd and ycor <= item 3 bd ] [
      set sext se0
      set blk i
    ]
  ]

  ask cells with [ blk >= 0 ] [
    set size 1.4
    set color blue
    stamp
    set size 0
  ]
end

; Active hyphae in star shape at centre of inoculum block (with tip at each end)
to setup-inoculum
  ask min-one-of cells [ distancexy ((cen-x0 + cen-x1) / 2) ((cen-y0 + cen-y1) / 2) ] [
    set si si0

    foreach dirs [ h ->
      let nb nbr-at h
      if nb != nobody [
        make-active nb h
        ask nb [ set si si0 ]
        make-tip nb h
      ]
    ]
  ]
end

to go
  do-step
  tick
  if simtime >= t-end or ticks > max-ticks [ stop ]
end

to do-step
  compute-dt
  move-tips
  branch
  uptake
  translocate
  maintain
  degrade
  diffuse-external

  ask cells [
    set si max (list 0 si)
    set sext max (list 0 sext)
  ]
  set simtime simtime + dt
end

; Adaptive timesteps so fast-moving tips don't skip over cells
to compute-dt
  let mr 0
  if any? tips [ set mr max [ tip-rate ] of tips ]
  let br b * max [ si ] of cells
  let dtv target-pmove / max (list mr br 1.0E-30)
  let diff max (list Di De Da-rich Da-poor)
  set dt min (list dtv (0.2 * delta-x * delta-x / diff))
end

to-report tip-rate
  let s [si] of tcell
  report 3 * Dp * s / (delta-x * delta-x) + v * s / delta-x
end

; Biased random walk: tips diffuse into three forward cells and convect straight ahead
; (With probabilities proportional to internal resource (Eqs. 1-2))
to move-tips
  ask tips [
    let s max (list 0 ([si] of tcell))
    let p-diff Dp * s * dt / (delta-x * delta-x)
    let p-conv v * s * dt / delta-x
    let move-total p-conv + 3 * p-diff

    ; Three forward directions: straight ahead (diffusion + convection), left and right (diffusion only)
    let routes (list (list heading (p-diff + p-conv))
                     (list (heading - 60) p-diff)
                     (list (heading + 60) p-diff))

    ; Drop directions blocked by the lattice boundary
    let open filter [ r -> (cell-toward (first r)) != nobody ] routes
    ifelse empty? open [
      die
    ][
      ; Redistribute blocked probability among open routes
      let share (move-total - sum map [ r -> last r ] open) / length open
      set open map [ r -> (list (first r) (last r + share)) ] open
      let stay 1 - move-total

      if random-float 1 >= stay [
        let h first rnd:weighted-one-of-list open [ r -> last r ]
        let nb cell-toward h

        ; Anastomosis: tip fuses if target cell already has hypha
        let fuse? any? [ my-links ] of nb

        ask tcell [
          make-active nb h
          set si si - c2 * delta-x
        ]

        ifelse fuse? [
          die
        ][
          set tcell nb
          set heading h
        ]
      ]
    ]
  ]
end

; New tips sprout from occupied cells with probability proportional to internal resource (rate b)
to branch
  refresh-tip-field
  ask cells with [ any? my-ahyphae and si > 0 and not tip-here? ] [
    if random-float 1 < b * si * dt [ branch-here ]
  ]
end

; Pick free direction beside existing hypha and branch
to branch-here
  ; Candidate directions: acute angles to growth direction of each active hypha
  let cand remove-duplicates reduce sentence map [ g ->
    (list ((g + 60) mod 360) ((g - 60) mod 360))
  ] growth-dirs

  ; Excluding directions already occupied
  let taken hypha-dirs
  set cand filter [ h -> not member? h taken and (nbr-at h) != nobody ] cand

  if empty? cand [ stop ]
  let d one-of cand
  let nb nbr-at d
  let tgt-occ any? [ my-links ] of nb

  make-active nb d
  set si si - c2 * delta-x
  if not tgt-occ [ make-tip nb d ]
end

; Hyphae absorb external resource (converting it to internal) (Eqs. 3-4)
to uptake
  ask cells [
    let m-hat (count my-ahyphae) * (delta-x / 2)
    let flux (max (list 0 si)) * m-hat * sext * dt
    set si si + c1 * flux
    set sext max (list 0 (sext - c3 * flux))
  ]
end

; Internal resource moves along active hyphae by diffusion and active transport towards tips (Eqs. 5-6)
to translocate
  refresh-tip-field
  ask cells [ set delta-si 0 ]

  ask ahyphae [
    let k end1
    let q end2
    let sk [si] of k
    let sq [si] of q

    ; Passive diffusion along the hypha
    let m-diff Di * dt * (sq - sk) / (delta-x * delta-x)

    ; Active transport towards tips (faster when resource-poor)
    let poor? ([sext] of k <= 0 and [sext] of q <= 0)
    let Da ifelse-value poor? [ Da-poor ] [ Da-rich ]
    let pgrad (ifelse-value [tip-here?] of k [1] [0]) - (ifelse-value [tip-here?] of q [1] [0])
    let m-act 0

    if pgrad > 0 [ set m-act Da * pgrad / delta-x * sq * dt / delta-x ]
    if pgrad < 0 [ set m-act Da * pgrad / delta-x * sk * dt / delta-x ]

    ask k [ set delta-si delta-si + m-diff + m-act ]
    ask q [ set delta-si delta-si - m-diff - m-act ]

    ; Metabolic cost of active transport
    let cost c4 * abs m-act
    if pgrad > 0 [ ask q [ set delta-si delta-si - cost ] ]
    if pgrad < 0 [ ask k [ set delta-si delta-si - cost ] ]
  ]

  ask cells [ set si si + delta-si ]
end

; Active hyphae consume internal resource
; Starved hyphae (si < omega) become inactive
to maintain
  ask cells [ set si si - (count my-ahyphae) * c5 * delta-x * dt ]
  ask cells with [ si < omega and any? my-ahyphae ] [
    ask my-ahyphae [ inactivate ]
    set si omega
  ]
end

; Inactive hyphae decay stochastically (Poisson process, rate d-i)
to degrade
  let prob (1 - exp (- d-i * dt))
  ask ihyphae [ if random-float 1 < prob [ die ] ]
end

; External resource diffuses within each resource block (Eq. 8)
; Only three of the six directions so each neighbouring pair is visited once
to diffuse-external
  ask cells [ set delta-se 0 ]

  ask cells with [ blk >= 0 ] [
    foreach (sublist dirs 0 3) [ h ->
      let nb nbr-at h
      if nb != nobody and [blk] of nb = blk [
        let m De * dt * (([sext] of nb) - sext) / (delta-x * delta-x)
        set delta-se delta-se + m
        ask nb [ set delta-se delta-se - m ]
      ]
    ]
  ]

  ask cells [ set sext sext + delta-se ]
end

; Mark which cells currently hold a growing tip
to refresh-tip-field
  ask cells [ set tip-here? false ]
  ask tips  [ ask tcell [ set tip-here? true ] ]
end

; The neighbouring cell lying in direction h (nobody at lattice edge)
to-report nbr-at [ h ]
  let hx xcor + sin h
  let hy ycor + cos h
  report one-of neighbours with [ distancexy hx hy < 0.25 ]
end

; The cell 1 lattice step ahead of a tip in direction h
to-report cell-toward [ h ]
  report [ nbr-at h ] of tcell
end

; Direction from this cell to another, snapped to the six lattice directions
to-report dir-to [ c ]
  report ((60 * round ((towards c - 30) / 60) + 30) mod 360)
end

; Directions of the hyphae leaving this cell
to-report hypha-dirs
  report map [ c -> dir-to c ] (sort link-neighbors)
end

; Growth directions of the active hyphae at this cell
to-report growth-dirs
  report map [ l -> [gdir] of l ] (sort my-ahyphae)
end

; Record a hypha between two cells (red for active, pink for inactive)
to make-active [ nb g ]
  if ihypha-neighbor? nb [ ask ihypha-with nb [ die ] ]
  if not ahypha-neighbor? nb [
    create-ahypha-with nb [ set color red set thickness 0.07 set gdir g mod 360 ]
  ]
end

to inactivate
  let nb end2
  ask end1 [ create-ihypha-with nb [ set color pink set thickness 0.06 ] ]
  die
end

to make-tip [ c h ]
  ask c [
    hatch-tips 1 [
      set tcell myself
      set heading h
      hide-turtle
    ]
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
210
10
658
459
-1
-1
8.0
1
10
1
1
1
0
1
1
1
0
54
0
54
0
0
1
ticks
30.0

BUTTON
21
41
88
74
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
118
41
181
74
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

MONITOR
21
119
95
164
t (sim time)
precision simtime 4
4
1
11

MONITOR
126
119
183
164
tips
count tips
17
1
11

MONITOR
9
185
96
230
active hyphae
count ahyphae
17
1
11

MONITOR
106
184
197
229
inactive hyphae
count ihyphae
17
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
