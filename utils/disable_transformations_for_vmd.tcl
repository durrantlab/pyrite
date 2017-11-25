set m {
 {{1.0 0.0 0.0 0.0}
  {0.0 1.0 0.0 0.0}
  {0.0 0.0 1.0 0.0}
  {0.0 0.0 0.0 1.0}}
}

set n [molinfo num]

for {set i 0} {$i < $n} {incr i} {
 molinfo ${i} set {center_matrix} $m
 molinfo ${i} set {rotate_matrix} $m
 molinfo ${i} set {scale_matrix} $m
 molinfo ${i} set {global_matrix} $m
}

