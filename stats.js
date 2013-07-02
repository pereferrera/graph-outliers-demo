function mainStats(graph, features, rankInfo, centrality, showCentralNodes, showCentralEdges, showInfluencers) {
	// For inspiration: http://www.slideshare.net/persuasion/facebook-network-analysis-using-gephi-11390089
    
    // Other metrics: http://en.wikipedia.org/wiki/Distance_(graph_theory)
		
	// Main stats
	
	var html = "";
	html += "<table class='table'><tbody>";
	html += "<tr><td># of nodes: </td><td>" + features.nNodes + "</td></tr>";
	html += "<tr><td># of edges: </td><td>" + features.edges + "</td></tr>";
	html += "<tr><td>Average degree: </td><td>" + features.averageDegree + "</td></tr>";
	html += "<tr><td>Max degree: </td><td>" + features.maxDegree + "</td></tr>";
	html += "<tr><td><a href='eigenvalue.html'>Largest eigenvalue</a>: </td><td>" + rankInfo.eigenvalue + "</td></tr>";
	html += "<tr><td># of disconnected communities: </td><td>" + findClusters(graph) + "</td></tr>";
	if(showCentralNodes) {
		html += "<tr><td>Diameter (longest shortest path): <td>" + centrality.diameter + "</td></tr>";
	}
	html += "<tr><td><a href='clustering.html'>Clustering coefficient</a>: <td>" + features.clusteringCoefficient + "</td></tr>";
	html += "</tbody>";
	
	$("#main-stats").html(html);
	
	// Influencers
	if(showInfluencers) {
		html = "<a href='influence.html'>Influencers:</a><ul>";
		for(var i = 0; i < 5 && i < rankInfo.ranks.length; i++) {
			html += "<li>" + rankInfo.ranks[i].node + " (" + rankInfo.ranks[i].rank + " e.centrality).</li>";
		}
		html += "</ul>";
		
		$("#influencers").html(html);
	}
	
	// Central nodes
	if(showCentralNodes) {
		centrality.centralities.sort(function(a, b) {
			return b.centrality - a.centrality;
		});
		
		html = "<a href='betweenness.html'>Central nodes:</a><ul>";
		for(var i = 0; i < 5 && i < centrality.centralities.length; i++) {
			html += "<li>" + centrality.centralities[i].node + " (" + centrality.centralities[i].centrality + " betweenness)</li>";
		}
		html += "</ul>";
	
		$("#central-nodes").html(html);
	}

	// Central edges
	if(showCentralEdges) {
		centrality.edgeCentralities.sort(function(a, b) {
			return b.centrality - a.centrality;
		});
		
		html = "<a href='betweenness.html'>Central edges:</a><ul>";
		for(var i = 0; i < 5 && i < centrality.edgeCentralities.length; i++) {
			html += "<li>" + centrality.edgeCentralities[i].from + " -> " + centrality.edgeCentralities[i].to + " (" + centrality.edgeCentralities[i].centrality + " betweenness)</li>";
		}
		html += "</ul>";
	
		$("#central-edges").html(html);
	}
}