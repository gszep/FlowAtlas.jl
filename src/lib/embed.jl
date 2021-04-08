function embed(data::DataFrame;path::AbstractString="",xdim::Int64=20,ydim::Int64=20,perplexity=10,maxIter=10000,eta=200.0)
	# Kratochvíl M, Koladiya A and Vondrášek J. Generalized EmbedSOM on quadtree-structured self-organizing maps. F1000Research 2020, 8:2120

	if isfile(path)
		som = deserialize(path)

	else
		som = initGigaSOM(data,xdim,ydim)
		som = trainGigaSOM(som,data)
		~isempty(path) && serialize(path,som)
	end

	########################### landmarks
	som.grid = tsne(som.codes,2,0,maxIter,perplexity;eta=eta)

	######################## extract clusters and embedding
	clusters = mapToGigaSOM(som,data)
	embedding = embedGigaSOM(som,data)

	rename!(clusters,"index"=>"cluster")
	return clusters, Matrix{Float64}(embedding')
end