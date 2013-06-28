function mainStats(graph, features, rankInfo, centrality) {
	// For inspiration: http://www.slideshare.net/persuasion/facebook-network-analysis-using-gephi-11390089
    
    // Other metrics: http://en.wikipedia.org/wiki/Distance_(graph_theory)
		
	// Main stats
	
	var html = "";
	html += "<table class='table'><tbody>";
	html += "<tr><td># of nodes: </td><td>" + features.nNodes + "</td></tr>";
	html += "<tr><td># of edges: </td><td>" + features.edges + "</td></tr>";
	html += "<tr><td>Average degree: </td><td>" + features.averageDegree + "</td></tr>";
	html += "<tr><td>Max degree: </td><td>" + features.maxDegree + "</td></tr>";
	html += "<tr><td>Largest eigenvalue: </td><td>" + rankInfo.eigenvalue + "</td></tr>";
	html += "<tr><td># of disconnected communities: </td><td>" + findClusters(graph) + "</td></tr>";
	html += "<tr><td>Diameter: <td>" + centrality.diameter + "</td></tr>";
	html += "<tr><td>Clustering coefficient: <td>" + features.clusteringCoefficient + "</td></tr>";
	html += "</tbody>";
	
	$("#main-stats").html(html);
	
	// Influencers
	
	html = "<ul>";
	for(var i = 0; i < 5; i++) {
		html += "<li>" + rankInfo.ranks[i].node + "</li>";
	}
	html += "</ul>";
	
	$("#influencers").html(html);

	// Central nodes
	
	centrality.centralities.sort(function(a, b) {
		return b.centrality - a.centrality;
	});
	
	html = "<ul>";
	for(var i = 0; i < 5; i++) {
		html += "<li>" + centrality.centralities[i].node + " (" + centrality.centralities[i].centrality + " betweenness)</li>";
	}
	html += "</ul>";

	$("#central-nodes").html(html);
	
	// Most-connected
	
	features.nodeFeatures.sort(function(a, b) {
		return b.degree - a.degree;
	});

	html = "<ul>";
	for(var i = 0; i < 5; i++) {
		html += "<li>" + features.nodeFeatures[i].node + " (" + features.nodeFeatures[i].degree + " mutual friends)</li>";
	}
	html += "</ul>";
}