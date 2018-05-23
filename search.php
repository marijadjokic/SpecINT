<!DOCTYPE html>
<html lang="en">
<head>
	<meta charset="utf-8">
	<title>SpecINT framework</title>
	<meta name="viewport" content="width=device-width, initial-scale=1.0">
	<meta name="description" content="Project Description">
	<meta name="author" content="Project Keywords">
	<link href="css/style.css" rel="stylesheet" type="text/css" />			
	<!--[if IE]><link href="css/style-ie.css" rel="stylesheet" type="text/css" /><![endif]-->	
	<script type="text/javascript" src="js/jquery-1.7.1.min.js"></script>
	<script type="text/javascript" src="js/site.js"></script>
	<!--[if IE]><script src="http://html5shiv.googlecode.com/svn/trunk/html5.js"></script><![endif]-->
	<script src="lib/codemirror.js"></script>
	<link rel="stylesheet" href="lib/codemirror.css">
	<script src="mode/sparql/sparql.js"></script>

	<script type="text/javascript">          
		$().ready(function () {

                // Setup form validation on the #register-form element
                $("#registration").validate({
                    // Specify the validation rules
                    rules: {
                        inchi: "required"
                    },
                    messages: {
                        inchi: "<br/>Please enter InChI"
                    }
                });
            });
       </script>
	   <style>
		   .letter:after {
				display: none;
			}
			
			.letter:before {
				display: none;
			}		
	   </style>
</head>
<body>

<div id="wrapper" class="letter">
	<section id="banner">
		<!-- <h1>How Do You <span>Nurture Your Money</span></h1> -->
		<img src="images/banner.png" alt="Banner">
	</section><!-- // end #banner -->
	<header id="header">
		<ul id="nav">
			<li class="active"><a href="index.html">HOME</a></li>
			<li><a href="example.html">EXAMPLES</a></li>
			<li><a href="search.php">ENDPOINT</a></li>
			<li><a href="references.html">REFERENCES</a></li>
			<li><a href="contact.html">SUPPORT</a></li>
		</ul>
	</header><!-- // end #header -->
	<div id="main-body" class="clearfix">
		<div id="content">
			<article> <!-- test -->
				<h1>SpecINT - user guide and run</h1>
				<p>				
					<b>SpecINT</b> provides to researchers the additional information about the same experimental results conducted in any independent
					laboratory without any prior knowledge about repository. User can use InChiKey notation for exploring data from online and centralized warehouse loaded into memory as a intermediate step for details listing.
				</p>
				<p>
					How to use software:
					<ol type="1">
						<li>Enter InChiKey for substance you want to use in experiments.</li>
						<li>Select external data source (e.g. CHEMBL) - the only data source at the moment.</li>
						<li>Select research topic.</li>
					</ol>
				</p>
				
				<div style="margin-bottom:20px;">
					For instance you can try following substances:<br /><br />
					<table class='results'>
					    <thead>
							<tr>
								<th scope='col'>Substance</th>
								<th scope='col'>InChIKey</th>
							</tr>
						</thead>
						<tbody>
							<tr>
								<td>Alfuzosin</td>
								<td>WNMJYKCGWZFFKR-UHFFFAOYSA-N</td>
							</tr>
							<tr>
								<td>Metixene</td>
								<td>MJFJKKXQDNNUJF-UHFFFAOYSA-N</td>
							</tr>
							<tr>
								<td>Cladribine</td>
								<td>PTOAARAWEBMLNO-KVQBGUIXSA-N</td>
							</tr>
							<tr>
								<td>3-(2-Amino-4-thiazolyl)-4-hydroxy-2H-1-benzopyran-2-one</td>
								<td>SCKPMNNSPZIZIF-UHFFFAOYSA-N</td>
							</tr>
							<tr>
								<td>Prazosin</td>
								<td>IENZQIKPVFGBNW-UHFFFAOYSA-N</td>
							</tr>
						</tbody>
					</table>
				</div>
				
				<div style="margin-bottom:20px;">
				    <img src="images/download.png" style="height:20px;" />
				    <a href="InChIKeys for experiments.xlsx">Download more examples</a>
				</div>
				
				<form name="registration" id="registration" method="post" >
				<fieldset class="field" style="position:relative;">
				<legend style="padding:0px 3px;font-weight:bold;">Explore substance:</legend>
					<table>
						<tr>
							<td>InChIKey*</td>
							<td>Information</td>
						</tr>
						<tr>
							<td>
								<input type="text" value="<?php echo isset($_POST['inchikey']) ? $_POST['inchikey'] : ''; ?>" name="inchikey" size="40"/>    
							</td>
							<td>
								<select id="information" name="information">
									<option value="target" <?php
									if (isset($_POST['information']) && (($_POST['information']) == 'target')) {
										echo('selected="selected"');
									}
									?>>Targets</option>
									<option value="cellline" <?php
									if (isset($_POST['information']) && (($_POST['information']) == 'cellline')) {
										echo('selected="selected"');
									}
									?>>Cell Lines</option>
									<option value="ic50" <?php
									if (isset($_POST['information']) && (($_POST['information']) == 'ic50')) {
										echo('selected="selected"');
									}
									?>>IC50 values</option>
								</select>
							</td>
						</tr>
						<tr >
							<td colspan=2 style="padding-top:15px;">
								<input type="submit" value="Submit" name="searching"/>                                
							</td>
						</tr>
					</table>
					</fieldset>
                </form>
				
				<p style="margin-top:20px;"><b style="color:red;">Information:</b> Good Internet connection is essential for software functioning. Calculation and links checking could take 5-10 minutes.</p> 
				
				<!--<div style="text-align:center;">
					<img src="http://cactus.nci.nih.gov/chemical/structure/InChIKey=WNMJYKCGWZFFKR-UHFFFAOYSA-N/image?footer=WNMJYKCGWZFFKR-UHFFFAOYSA-N&width=500" />
				</div>-->
				
				<?php
                    set_time_limit(0);   
					include_once('./semsol-arc2/semsol/ARC2.php'); /* ARC2 static class inclusion */
					 
					$config = array(/* remote endpoint configuration */
						'remote_store_endpoint' => 'http://localhost:3030/PIBASTEST/sparql',
					);

					$store = ARC2::getRemoteStore($config); /* instantiation */
				
					if (isset($_POST['searching'])) {
						if ((isset($_POST['inchikey'])) && (isset($_POST['information']))) {
							$inchikey = $_POST['inchikey'];
							$information = $_POST['information'];
							$command = escapeshellcmd("python ./search.py " . $inchikey ." ".$information);						
							$output = shell_exec($command);
							//echo $output;
							
							$for_fv_and_pr = substr($output, 0, strpos($output, '}') + 1);
							$fv = substr($for_fv_and_pr, strpos($for_fv_and_pr, '['), strpos($for_fv_and_pr, ']') - 1);
							$pr = substr($for_fv_and_pr, strpos($for_fv_and_pr, '{'), strpos($for_fv_and_pr, '}'));
							$for_addtional_imformation = substr($output, strpos($output, '}') + 2);
							$for_addtional_imformation_new = substr($output, strpos($output, '}') + 2);
							$additional_information = explode(" ", $for_addtional_imformation);
							
							//var_dump($fv);

							if ((strpos($output, 'fv') !== false) and (strpos($output, 'pr')  === false) and (strpos($output, $information)  === false)) {
								$fv = substr($output, strpos($output, '-')+1);
								echo '<table>
										  <tr>
											  <td><img style="width: 300px ;height: 200px" src="img/completed_graph_'.$_POST['inchikey'].'.png"/></td>
										  </tr>
										  <tr>
										   <td>' . $fv . '</td> 
										  </tr>
									  </table>';
							}

							if ((strpos($output, 'fv') !== false) and (strpos($output, 'pr') !== false) and (strpos($output, $information)  === false)) {
								$fv = substr($output, strpos($output, '['), strpos($output, ']') - 1);
								$for_pr=substr($output, strpos($output, ']')+1);
								if(strpos($for_pr, '{')!=''){
								  $pr=substr($for_pr, strpos($for_pr, '{'), strpos($for_pr, '}')-2);
								}
								else
								{
									$pr= substr($for_pr, strpos($for_pr, '-')+1);
								}
								
								echo '<table>
										  <tr>
											  <td><img style="width: 300px ;height: 200px" src="img/completed_graph_'.$_POST['inchikey'].'.png"/></td>
											  <td><img style="width: 300px ;height: 200px" src="img/noncompleted_graph_'.$_POST['inchikey'].'.png"/></td>
										  </tr>
										  <tr>
											  <td>' . $fv . '</td>                
											  <td>' . $pr . '</td>
										  </tr>
									  </table>
									  ';
								echo '<table><tr><td>"No additional results!"</td></tr></table>';
							}

							if ((strpos($output, 'fv') !== false) and (strpos($output, 'pr') !== false) and (strpos($output, $information) !== false)) {
								$for_fv_and_pr = substr($output, 0, strpos($output, '}') + 1);
								$fv = substr($for_fv_and_pr, strpos($for_fv_and_pr, '['), strpos($for_fv_and_pr, ']') - 1);
								$pr = substr($for_fv_and_pr, strpos($for_fv_and_pr, '{'), strpos($for_fv_and_pr, '}'));
								
								$for_addtional_imformation1 = substr($output, strpos($output, '}') + 2);
                                #echo $for_addtional_imformation1;
                                $for_addtional_imformation2 = substr($for_addtional_imformation1, strpos($for_addtional_imformation1, '-') + 2, strpos($for_addtional_imformation1, '^') - 8);
                                #echo $for_addtional_imformation2;
                                $additional_information = explode(" ", $for_addtional_imformation2);
                                $query = substr($output, strpos($output, '^') + 1);
								
								/*$for_addtional_imformation = substr($output, strpos($output, '}') + 2);
								$for_addtional_imformation_new = substr($output, strpos($output, '}') + 2);
								$additional_information = explode(" ", $for_addtional_imformation);
								var_dump($additional_information);*/

								/*echo '<table>
										  <tr>
											  <td><img style="width: 300px ;height: 200px" src="img/completed_graph_'.$_POST['inchikey'].'.png"/></td>
											  <td><img style="width: 300px ;height: 200px" src="img/noncompleted_graph_'.$_POST['inchikey'].'.png"/></td>
										  </tr>
										  <tr>
											  <td>' . $fv . '</td>                
											  <td>' . $pr . '</td>
										  </tr>
									  </table>
									  ';*/
								  
							   if($additional_information!=" "){
									echo "<textarea id='queryString' rows='10' cols='100'>"; 
										echo $query;
									echo "</textarea>";
									
									// crtanje supstance
									echo "<div style='text-align:center;margin-bottom:20px;'>";
									echo "<img src='http://cactus.nci.nih.gov/chemical/structure/InChIKey=".$inchikey."/image?footer=".$inchikey."&width=500' />";
									echo "</div>";

									echo "<table class='results'>";
									echo "<thead>";
										
									if ($information == 'target') {
										echo "<tr>
												<th scope='col'>Target</th>
												<th scope='col'>Target Report Card</th>
											 </tr>";
										// zatvori head i otvori body
										echo "</thead>
											<tbody>";
											
										for ($x = 0; $x <= count($additional_information) - 2; $x++) {

											if (filter_var($additional_information[$x], FILTER_VALIDATE_URL)) {
												if (strpos($additional_information[$x], 'chembl/target') !== false) {
													$link_array = explode('/', $additional_information[$x]);
													$id = end($link_array);
													echo '<tr><td><a href=' . $additional_information[$x] . ' target="_blank">' . $additional_information[$x] . '</a></td>';
													echo '<td><b style="margin-right:10px;">target card</b> <a href="https://www.ebi.ac.uk/chembldb/index.php/target/inspect/' . $id . '" target="_blank"><img src="images/report.png" width="50" style="margin-right:10px;" /></a></td></tr>';
												}

												if (strpos($additional_information[$x], 'bio2rdf.org/drugbank') !== false) {
													$link_array = explode(':', $additional_information[$x]);
													$id = end($link_array);
													echo '<tr><td><a href=' . $additional_information[$x] . ' target="_blank">' . $additional_information[$x] . '</a></td>';
													echo '<td><b style="margin-right:10px;">target card</b> <a href="http://www.drugbank.ca/biodb/bio_entities/' . $id . '" target="_blank"><img src="images/report.png" width="50" style="margin-right:10px;" /></a></td></tr>';
												}

												if (strpos($additional_information[$x], 'chemogenomics/resource') !== false) {
													$link_array = explode(':', $additional_information[$x]);
													$id = end($link_array);
													echo '<tr><td><a href=' . $additional_information[$x] . ' target="_blank">' . $additional_information[$x] . '</a></td>';
													echo '<td><b style="margin-right:10px;">target card</b> <a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=' . $id . '" target="_blank"><img src="images/report.png" width="50" style="margin-right:10px;" /></a></td></tr>';
												}
											} else {
												echo '<tr><td>' . $additional_information[$x] . '</td><td>no more information</td></tr>';
											}
										}
									}
										
									if ($information == 'cellline') {
										echo "<tr>
												<th scope='col'>Cell line</th>
												<th scope='col'>Cell Line Report Card</th>
											 </tr>";
										// zatvori head i otvori body
										echo "</thead>
											<tbody>";	 
										
										for ($x = 0; $x <= count($additional_information) - 2; $x++) {
											if (filter_var($additional_information[$x], FILTER_VALIDATE_URL)) {
												if (strpos($additional_information[$x], 'chembl/cell_line') !== false) {
													$link_array = explode('/', $additional_information[$x]);
													$id = end($link_array);
													echo '<tr><td><a href=' . $additional_information[$x] . ' target="_blank">' . $additional_information[$x] . '</a></td>';
													echo '<td><b style="margin-right:10px;">cell-line card</b> <a href="https://www.ebi.ac.uk/chembldb/cell/inspect/' . $id . '" target="_blank"><img src="images/report.png" width="50" style="margin-right:10px;" /></a></td></tr>';
												}
											} 
											else {
												echo '<tr><td>' . $additional_information[$x] . '</td><td>no more information</td></tr>';
											}
										}
									}	
									
									if ($information == 'ic50') {
										echo "<tr>
												<th scope='col'>Cell line</th>
												<th scope='col'>IC<sub>50</sub> value</th>
											 </tr>";
										// zatvori head i otvori body
										echo "</thead>
											<tbody>";
											
										for ($x = 0; $x <= count($additional_information) - 2; $x++) {
											if ($x % 2 == 0)
											{								
											    if (filter_var($additional_information[$x], FILTER_VALIDATE_URL)) {											    
													if (strpos($additional_information[$x], 'chembl/cell_line') !== false){
														$link_array = explode('/', $additional_information[$x]);
														$id = end($link_array);											
														echo '<tr><td><a href="https://www.ebi.ac.uk/chembl/cell/inspect/'.$id.'" target="_blank"> https://www.ebi.ac.uk/chembl/compound/inspect/'.$id.'</a></td>';											   
													}
												}
												else
												{
													echo '<tr><td>' . $additional_information[$x] . '</td>';											   
												} 
											}
											else
											{									  
												echo '<td>' . $additional_information[$x] . '</td></tr>';								  
											}										
										}
									}

									echo "</tbody></table>";
									//echo '<b>Query:</b><br/>' . nl2br(htmlspecialchars($query));
								}else {
									echo 'No additional information!';
								}
								
								echo "<h1>Graph background</h1>";
								echo "Here you can see potential edges between data sources and existing edges at the moment.";
								echo '<table style="margin-top:30px;">
										<tr>
											<td><img style="width:400px;" src="img/completed_graph_'.$_POST['inchikey'].'.png"/></td>
											<td><b>Fiedler eigenvector:</b><br /><br />' . $fv . '</td>
											
										</tr>
										<tr>
										    <td><img style="width:400px;" src="img/noncompleted_graph_'.$_POST['inchikey'].'.png"/></td>                
										    <td><b>PageRank eigenvector:</b><br /><br />' . $pr . '</td>
										</tr>
									</table>';
							}                          

							if ((strpos($output, 'fv-') === false) and (strpos($output, 'pr-')  === false) and (strpos($output, $information)  === false)) {
								echo $output;
							}
						} else {
                            echo '<script type="text/javascript">alert("Please, enter/select data.")</script>';
                        }						
					}
					?>
			</article>
		</div><!-- // end #content -->
		<aside id="sidebar">
			<h1>LATEST NEWS</h1>
            <a class="twitter-timeline"  href="https://twitter.com/spectralni/lists/bioinformatics1" data-widget-id="658566232925409280">Tweets from https://twitter.com/spectralni/lists/bioinformatics1</a>
            <script>!function(d,s,id){var js,fjs=d.getElementsByTagName(s)[0],p=/^http:/.test(d.location)?'http':'https';if(!d.getElementById(id)){js=d.createElement(s);js.id=id;js.src=p+"://platform.twitter.com/widgets.js";fjs.parentNode.insertBefore(js,fjs);}}(document,"script","twitter-wjs");</script>                          
		</aside><!-- // end #sidebar -->
	</div><!-- // end #main-body -->
	<footer id="footer">
		<!-- <span class="fl"><a href="#">Terms of Use</a> | <a href="#">Privacy Policy</a></span> -->
		<span class="fr">Copyright &copy; <strong>SpecINT</strong>. All right reserved.</span>
	</footer><!-- // end #footer -->
</div><!-- // end #wrapper -->
<p class="author"></p>
<script>
	  var editor = CodeMirror.fromTextArea(document.getElementById("queryString"), {
	  mode: "application/sparql-query",
	  matchBrackets: true,
	  lineNumbers: true,
	  readOnly: true
	  });
</script>
</body>
</html>