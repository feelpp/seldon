<?php echo '<?xml version="1.0"  encoding="iso-8859-1"?'.'>' ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">

<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

<?php $root='..';?>

<head>
<title>Seldon user's guide</title>
<link rel="stylesheet" type="text/css" href="<?php echo $root?>/content.css"/>
<link rel="stylesheet" href="tabs.css" type="text/css"/>
<link rel="stylesheet" href="guide.css" type="text/css"/>
<script type="text/javascript" src="<?php echo $root.'/';?>prettify.js"></script>
</head>

<body onload="prettyPrint()">

<div class="page">

<?php include $root.'/header.php';?>

<div class="doc">

<?php function HL($file_, $section_, $string_)
{
if ($file_ == $section_)
  echo '<em>'.$string_.' </em>';
else
  echo '<a href="'.$section_.'.php">'.$string_.'</a>';
}; ?>

<?php $file=basename($file, ".php"); $file = explode(".", $file); $file = $file[0];?>

<div class="nav">

<ul>
<li class="jelly"> <b>USER'S GUIDE</b> </li>
<li class="jelly"> <?php HL($file, "installation", "Installation");?> </li>
<li class="jelly"> <?php HL($file, "overview", "Overview");?> </li>
<li class="jelly"> <?php HL($file, "vectors", "Vectors");?>

<?php if (basename($_SERVER['REQUEST_URI'], ".php") == "vectors"
or basename($_SERVER['REQUEST_URI'], ".php") == "class_vector"
or basename($_SERVER['REQUEST_URI'], ".php") == "class_sparse_vector"
or basename($_SERVER['REQUEST_URI'], ".php") == "functions_vector")
{
  echo '<ul class="navsubul"> <li class="jelly">';
  HL($file, "class_vector", "Dense Vectors");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "class_sparse_vector", "Sparse Vectors");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "functions_vector", "Functions");
  echo '</li> </ul>';
} ?>

</li>

<li class="jelly"> <?php HL($file, "matrices", "Matrices");?>

<?php if (basename($_SERVER['REQUEST_URI'], ".php") == "matrices"
or basename($_SERVER['REQUEST_URI'], ".php") == "class_matrix"
or basename($_SERVER['REQUEST_URI'], ".php") == "class_sparse_matrix"
or basename($_SERVER['REQUEST_URI'], ".php") == "functions_matrix")
{
  echo '<ul class="navsubul"> <li class="jelly">';
  HL($file, "class_matrix", "Dense Matrices");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "class_sparse_matrix", "Sparse Matrices");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "functions_matrix", "Functions");
  echo '</li> </ul>';
} ?>

</li>

<li class="jelly"> <?php HL($file, "array3d", "3D&nbsp;Arrays");?>  </li>
<li class="jelly"> <?php HL($file, "allocators", "Allocators");?>  </li>
<li class="jelly"> <?php HL($file, "exceptions", "Exceptions");?>  </li>
<li class="jelly"> <?php HL($file, "computations", "Computations");?>

<?php if (basename($_SERVER['REQUEST_URI'], ".php") == "computations"
or basename($_SERVER['REQUEST_URI'], ".php") == "functions_blas"
or basename($_SERVER['REQUEST_URI'], ".php") == "functions_lapack"
or basename($_SERVER['REQUEST_URI'], ".php") == "direct"
or basename($_SERVER['REQUEST_URI'], ".php") == "iterative")
{
  echo '<ul class="navsubul"> <li class="jelly">';
  HL($file, "functions_blas", "Blas");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "functions_lapack", "Lapack");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "direct", "Direct Solvers");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "iterative", "Iterative Solvers");
  echo '</li> </ul>';
} ?>

</li>

<li class="jelly"> <?php HL($file, "python", "Python Interface");?> </li>
<li class="jelly"> <?php HL($file, "glossary", "Index");?> </li>
<li class="jelly"> <b>API REFERENCE</b> </li>
<li class="jelly"> <?php HL($file, "annotated", "Classes");?>
<ul class="navsubul"> <li class="jelly"> <?php HL($file, "annotated", "Class List");?> </li> 
<li class="jelly"> <?php HL($file, "hierarchy", "Class Hierarchy");?> </li>
<li class="jelly"> <?php HL($file, "functions", "Class Members");?>
</li> </ul> </li>
<li class="jelly"> <?php HL($file, "namespacemembers", "Functions");?> </li>
<li class="jelly"> Search for <form action="search.php" method="get">
    <input class="search" type="text" name="query" value="" size="20" accesskey="s"/>
  </form>
</li>
<!-- <li class="jelly"> <?php HL($file, "faq", "F.A.Q.");?> </li>-->
<li class="jelly"> <a
href="mailto:Marc.Durufle@math.u-bordeaux1.fr,Vivien.Mallet@inria.fr"
style="color:black">Support</a></li>
</ul>

</div>

<div class="doxygen">
