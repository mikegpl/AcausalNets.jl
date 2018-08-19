#=
sets:
- Julia version: 
- Author: marcin
- Date: 2018-08-19
=#

is_subset(s1::Set, s2::Set)::Bool = intersect(s1, s2) == s1