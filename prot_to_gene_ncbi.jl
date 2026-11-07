#!/usr/bin/env julia -t auto --project=@.

using CodecZlib
using CSV
using DataFrames
using GFF3
using TranscodingStreams

@assert length(ARGS) >= 2 "Usage: julia script.jl <input_gff3.gz> <output_csv>"

input_path = ARGS[1]
output_path = ARGS[2]

if endswith(input_path, ".gz")
    strm = TranscodingStream(GzipDecompressor(), open(ARGS[1]))
    rdr = GFF3.Reader(strm)
else
    rdr = open(GFF3.Reader, input_path)
end

recs = [rec for rec in rdr]

close(rdr)

all_ids = []

for rec in recs

    atr_dict = GFF3.attributes(rec) |> Dict

    if !haskey(atr_dict, "Dbxref")
        continue
    end
    
    idx = findfirst(atr -> startswith(atr, "GeneID:"), atr_dict["Dbxref"])
    
    if isnothing(idx)
        continue
    end
    
    push!(all_ids, atr_dict["Dbxref"][idx])
end

all_ids = filter(!isnothing, all_ids) |> unique
id_dict = [id => [] for id in all_ids] |> Dict

for rec in recs
    atr_dict = GFF3.attributes(rec) |> Dict
    
    if !haskey(atr_dict, "Dbxref") || !haskey(atr_dict, "ID")
        continue
    end
    
    idx = findfirst(atr -> startswith(atr, "GeneID:"), atr_dict["Dbxref"])
    
    if isnothing(idx)
        continue
    end
    
    parent_id = atr_dict["Dbxref"][idx]
    
    push!(id_dict[parent_id], first(atr_dict["ID"]))
end

prot_to_gene = Dict{String, String}()
for child_ids in values(id_dict)
    gene_idx = findfirst(atr -> startswith(atr, "gene-"), child_ids)
    if isnothing(gene_idx)
        continue
    end

    prot_idxs = findall(atr -> startswith(atr, "cds-"), child_ids)
    if isempty(prot_idxs)
        continue
    end

    for prot_idx in prot_idxs
        prot_id = child_ids[prot_idx]
        gene_id = child_ids[gene_idx]
        prot_to_gene[prot_id] = gene_id
    end
end

ps = [replace(p[1], "cds-" => "") for p in prot_to_gene |> collect]
gs = [replace(p[2], "gene-" => "") for p in prot_to_gene |> collect]
df = DataFrame(:gene => gs, :protein => ps)
CSV.write(output_path, df)
