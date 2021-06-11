function embed(data::DataFrame;path::AbstractString="",xdim::Int64=20,ydim::Int64=20,perplexity=10,maxIter=10000,eta=200.0)
	# Kratochvíl M, Koladiya A and Vondrášek J. Generalized EmbedSOM on quadtree-structured self-organizing maps. F1000Research 2020, 8:2120

	if isfile(path)
		@info """Using embedding saved at $path"""
		som = deserialize(path)

	else
		@info """Calculating embedding..."""
		som = initGigaSOM(data,xdim,ydim)
		som = trainGigaSOM(som,data)
		~isempty(path) && serialize(path,som)
	end

	######################## extract clusters and embedding
	try 
		som.grid = tsne(som.codes,2,0,maxIter,perplexity;eta=eta)
		clusters = mapToGigaSOM(som,data)
		embedding = embedGigaSOM(som,data)

		rename!(clusters,"index"=>"cluster")
		return clusters, Matrix{Float64}(embedding')
	catch

		@info """Re-calculating embedding..."""
		som = initGigaSOM(data,xdim,ydim)
		som = trainGigaSOM(som,data)
		~isempty(path) && serialize(path,som)

		som.grid = tsne(som.codes,2,0,maxIter,perplexity;eta=eta)
		clusters = mapToGigaSOM(som,data)
		embedding = embedGigaSOM(som,data)

		rename!(clusters,"index"=>"cluster")
		return clusters, Matrix{Float64}(embedding')
	end
end