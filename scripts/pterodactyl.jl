"""
# Pterodactyl
This module contains a list of functions
useful for working with NAUTILUS network files
## Main Functions in Pterodactyl
* `Pterodactyl.searchNetwork()`
* `Pterodactyl.printReactions()`
* `Pterodactyl.makenew()`
* `Pterodactyl.writeReactions()`
## Help
Type `?` followed by the function name for help. For example:
`?Pterodactyl.searchNetwork`
# Typical function order
I recommend calling functions in the following order
1. `Pterodactyl.searchNetwork()`
2. `Pterodactyl.makenew()`
"""
module 
Pterodactyl
using DataFrames
using Formatting

path = pwd()

println("Welcome to Pterodactyl! Type \"?Pterodactyl\" for help.")

"Get input from user"
function input(prompt::AbstractString="")
  println(prompt)
  return chomp(readline())
end

"Take NAUTILUS style species name and add HTML-style super/sub-scripts"
function formatMolecule(specieslist)
  i = 1
  for species in specieslist
    species = replace(species,"+","<sup>+</sup>")
    species = replace(species,"-","<sup>-</sup>")
    species = replace(species,"*","<sup>*</sup>") 
    for letter in species
      test = tryparse(Int64,"$letter")
      if ( isnull(test) == false )
        species = replace(species,"$letter","<sub>$letter</sub>")
      end
    end
    specieslist[i] = species
    i += 1
  end
  return specieslist
end

"""
# searchNetwork
Search a \".in\" network file for all reactions involving a certain species
# Output
* `prod_df`: Production reactions
* `dest_df`: Destruction reactions
* `both_df`: Both production and destruction reactions
"""
function searchNetwork()
  # Select gas/grain reactions
  global gasgrain = ""
  while true
    gasgrain = input("Enter (a) for gas species or (r) for grain species : ")
    if gasgrain in ["a","r"]
      break 
    else
      println("Try again: ")
    end
  end

  if gasgrain == "a"
    species_file = "gas_species.in"
    global reactions_file = "gas_reactions.in"
    gasgrain = "gas"
  else
    species_file = "grain_species.in"
    global reactions_file = "grain_reactions.in"
    gasgrain = "grain"
  end

  # Check whether species is in list
  toopen = path*"/"*species_file
  println("Now opening $species_file")
  species_df = readtable(
                         path*"/"*species_file,
                         separator = ' ',
                         names = [
                                  :species, 
                                  :charge, 
                                  :H, 
                                  :He, 
                                  :C, 
                                  :N, 
                                  :O, 
                                  :Si, 
                                  :S, 
                                  :Fe, 
                                  :Na,
                                  :Mg, 
                                  :Cl, 
                                  :P, 
                                  :F
                                 ],
                         eltypes = [
                                    String,
                                    Int64,
                                    Int64,
                                    Int64,
                                    Int64,
                                    Int64,
                                    Int64,
                                    Int64,
                                    Int64,
                                    Int64,
                                    Int64,
                                    Int64,
                                    Int64,
                                    Int64,
                                    Int64
                                   ],
                         header = false
                        )


  # Get species name
  global species = ""
  while true
    species = input("Enter species name: ")
    if species in species_df[:species] 
      println("You chose $species")
      break 
    else
      println("Try again: ")
    end
  end

  # Get reactions in DataFrame
  pdb = "b"
  #=
  while true
    pdb = input("Production, Destruction, or Both? (p/d/b) ")
    if pdb in ["p","b","d"]
      break 
    else
      println("Try again: ")
    end
  end
  =#

  # Make reactions dataframes
  global reactions_df = DataFrame(
                                  r1 = String[],
                                  r2 = String[],
                                  r3 = String[],
                                  p1 = String[],
                                  p2 = String[],
                                  p3 = String[],
                                  p4 = String[],
                                  p5 = String[],
                                  α  = Float64[],
                                  β  = Float64[],
                                  γ  = Float64[],
                                  fx1 = Float64[],
                                  fx2 = Float64[],
                                  sx3 = String[],
                                  itype = Int64[],
                                  tmin = Float64[],
                                  tmax = Float64[],
                                  formula = Int64[],
                                  id = Int64[],
                                  ix4 = Int64[],
                                  ix5 = Int64[]
                                 )

  global prod_df = DataFrame() 
  global dest_df = DataFrame() 

  # Load reactions file
  infile = open(pwd()*"/"*reactions_file)
  readline(infile)
  i = 0
  println("Now opening $reactions_file")
  for ln in eachline(infile)
    if ln[1:1] == "!"
      continue
    elseif ln[1:1] == " "
      continue
    end
    i += 1
    if i%1000 == 0
      println("On line $i")
    end
    try
      push!(
            reactions_df,
            [
             strip(ln[1:11]),  #r1
             strip(ln[12:22]), #r2
      strip(ln[23:34]), #r3
      strip(ln[35:45]), #p1
      strip(ln[46:56]), #p2
      strip(ln[57:67]), #p3
      strip(ln[68:78]), #p4
      strip(ln[79:90]), #p5
      parse(Float64,ln[91:101]), #α
      parse(Float64,ln[102:112]), #β
      parse(Float64,ln[113:123]), #γ
      parse(Float64,ln[124:132]), #fx1
      parse(Float64,ln[133:141]), #fx2
      strip(ln[142:146]), #sx3
      parse(Int64,ln[147:150]), #itype
      parse(Float64,ln[151:157]), #tmin
      parse(Float64,ln[158:163]), #tmax
      parse(Int64,ln[164:166]), #formula
      parse(Int64,ln[167:172]), #id
      parse(Int64,ln[173:175]), #ix4
      parse(Int64,ln[176:end]), #ix5
     ]
           )
    catch y
      if isa(y,BoundsError)
        return
      end
    end
  end

  # Check for destruction
  if pdb in ["d","b"]
    dest_df = vcat(
                   dest_df,
                   reactions_df[reactions_df[:r1].==species,:],
                   reactions_df[reactions_df[:r2].==species,:],
                  )
  end

  # Check for production
  if pdb in ["p","b"]
    prod_df = vcat(
                   prod_df,
                   reactions_df[reactions_df[:p1].==species,:],
                   reactions_df[reactions_df[:p2].==species,:],
                   reactions_df[reactions_df[:p3].==species,:],
                   reactions_df[reactions_df[:p4].==species,:],
                  )

  end


  global both_df = DataFrame()
  both_df = vcat(both_df, prod_df, dest_df)
  println()
  println("The production reactions are in `prod_df`")
  println("The destruction reactions are in `dest_df`")
  println("Both are in `both_df`")
end

"Nicely formatted printing of reactions"
function printReactions(out_df)
  #Print reactions
  fmt = "%-11c + %-11c -> %-11c + %-11c + %-11c"
  println()
  for i in 1:length(out_df)
    #    s = sprintf1(fmt,out_df[:r1][i],out_df[:r2][i],out_df[:p1][i],out_df[:p2][i],out_df[:p3][i])
    try
      x = out_df[:r1][i]
    catch y
      if isa(y,BoundsError)
        return
      end
    end
    if out_df[:r2][i] == ""
      if out_df[:p3][i] != ""
        println("$i: $(out_df[:r1][i]) -> $(out_df[:p1][i]) + $(out_df[:p2][i]) + $(out_df[:p3][i])")
      elseif out_df[:p2][i] != ""
        println("$i: $(out_df[:r1][i]) -> $(out_df[:p1][i]) + $(out_df[:p2][i])")
      else
        println("$i: $(out_df[:r1][i]) -> $(out_df[:p1][i])")
      end

    else
      if out_df[:p3][i] != ""
        println("$i: $(out_df[:r1][i]) + $(out_df[:r2][i]) -> $(out_df[:p1][i]) + $(out_df[:p2][i]) + $(out_df[:p3][i])")
      elseif out_df[:p2][i] != ""
        println("$i: $(out_df[:r1][i]) + $(out_df[:r2][i]) -> $(out_df[:p1][i]) + $(out_df[:p2][i])")
      else
        println("$i: $(out_df[:r1][i]) + $(out_df[:r2][i]) -> $(out_df[:p1][i])")
      end
    end
  end
end

"""
Write new reactions to gas/grain reactions file
# Arguments
* `out_df::DataFrames.DataFrame`: The reactions to write to file
"""
function writeReactions(out_df)
  # Select gas/grain reactions
  gasgrain = ""
  ggr = ""
  while true
    gasgrain = input("Enter (a) for gas species or (r) for grain species : ")
    if gasgrain in ["a","r"]
      if gasgrain == "a"
        ggr = "gas"
      else
        ggr = "grain"
      end
      while true
        answer = input("You selected $ggr\nis this correct? (y/n) : ")
        if answer in ["y","n"]
          if answer == "y"
            break
          end
        else
          println("Try again: ")
        end
      end
      break
    else
      println("Try again: ")
    end
  end

  reactions_file = "$(ggr)_reactions.in"

  # Get last reaction ID in file
  outfile = open(pwd()*"/"*reactions_file,"r+")
  last_id = 0
  for ln in eachline(outfile)
    if !(ln[1:1] in ["!"," "])
      last_id = parse(Int64,ln[167:172]) #id
    end
  end
  println("The last ID num is $last_id")

  write(outfile,"!\n")
  write(outfile,"!***** Reactions below added by pterodactyl.jl *****\n")
  for i in 1:length(out_df)
    try
      x = out_df[:r1][i]
    catch y
      if isa(y,BoundsError)
        return
      end
    end
    last_id += 1
    line = @sprintf("%-11s%-11s%-12s%-11s%-11s%-11s%-11s%-11s%11.3E%11.3E%11.3E%9.2E%9.2E%5s%3d%7d%7d%3d%6d%2d%3d\n",
                    out_df[:r1][i],
                    out_df[:r2][i],
                    out_df[:r3][i],
                    out_df[:p1][i],
                    out_df[:p2][i],
                    out_df[:p3][i],
                    out_df[:p4][i],
                    out_df[:p5][i],
                    out_df[:α][i],
                    out_df[:β][i],
                    out_df[:γ][i],
                    out_df[:fx1][i],
                    out_df[:fx2][i],
                    out_df[:sx3][i],
                    out_df[:itype][i],
                    out_df[:tmin][i],
                    out_df[:tmax][i],
                    out_df[:formula][i],
                    last_id,
                    out_df[:ix4][i],
                    out_df[:ix5][i],
                   )
    write(outfile,line)
  end
  close(outfile)
  #= 
  If you've just added/renumbered the gas-phase reactions, then you need to 
  renumber the grain-surface reactions as well. Open the grain-surface reactions
  file and renumber those.
  =# 
  if reactions_file == "gas_reactions.in"
    last_num = renumber(last_id,"grain_reactions.in")
  end
end

"Renumber reaction IDs as needed in \".in\" files"
function renumber(num,infilename)
  println("Now renumbering grain_reactions")
  infile = open(pwd()*"/"*infilename,"r")
  outfile = open(pwd()*"/renumbered_"*infilename,"w")
  num += 1
  for line in eachline(infile)
    if line[1:1] == "!"
      continue
    end
    global id = line[167:172]
    try
      # Try to convert the appropriate line section to an int
      parse(Int64,id) #id
    catch
      # If there's an error, just continue on
    end
    renum = @sprintf("%6d",num)
    newline = line[1:165]*renum*line[172:end] 
    println(line)
    println(newline)
    write(outfile,newline)
    num += 1
  end
  close(outfile)
  close(infile)
  println("At the end of this subroutine, last_num = $num")
  return num
end

"""
Replace one product/reactant in a set of reactions with another species
input by the user. The output of this function is the same set of reactions
with one species replaced.
# Arguments
* `reactions::DataFrames.DataFrame`: A set of reactions
"""
function makenew(reactions)
  # Get species to replace
  toreplace = ""
  while true
    toreplace = input("What am I looking for? : ")
    answer = input("You selected $toreplace \nis this correct? (y/n): ")
    if answer in ["y","n"] 
      if answer == "y"
        break 
      end
    else
      println("Try again: ")
    end
  end
  # Get new species
  replacewith = ""
  while true
    replacewith = input("What am I replacing it with? : ")
    answer = input("You selected $replacewith \nis this correct? (y/n): ")
    if answer in ["y","n"] 
      if answer == "y"
        break 
      end
    else
      println("Try again: ")
    end
  end

  out_df = DataFrame()
  out_df = reactions
  for i in 1:size(reactions,1)
    if reactions[:r1][i] == toreplace  
      out_df[:r1][i] = replacewith
    end
    if reactions[:r2][i] == toreplace
      out_df[:r2][i] = replacewith
    end
    if reactions[:r3][i] == toreplace
      out_df[:r3][i] = replacewith
    end
    if reactions[:p1][i] == toreplace
      out_df[:p1][i] = replacewith
    end
    if reactions[:p2][i] == toreplace
      out_df[:p2][i] = replacewith
    end
    if reactions[:p3][i] == toreplace
      out_df[:p3][i] = replacewith
    end
  end

  # Print reactions
  answer = ""
  while true
    answer = input("View new reactions? (y/n): ")
    if answer in ["y","n"] 
      if answer == "y"
        break 
      end
    else
      println("Try again: ")
    end
  end
  if answer == "y"
    printReactions(out_df)
  end


  # Write reactions to file
  writeyn = ""
  while true
    writeyn = input("Write reactions to file? (y/n) : ")
    answer = input("You selected $writeyn\nis this correct? (y/n): ")
    if answer in ["y","n"] 
      if answer == "y"
        break 
      end
    else
      println("Try again: ")
    end
  end

  if writeyn == "y"
    writeReactions(out_df)
  end

  global new_reactions = DataFrame()
  new_reactions = out_df
  println("Type `new_reactions` to view new reactions again,\n
          or print again using `printReactions`")
end
end
