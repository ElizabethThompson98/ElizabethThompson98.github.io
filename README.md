<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ElizabethThompson98</title>

    <style>
        body {
            margin: 20px; /* Add some margin to the body to provide space */
            display: flex;
            flex-wrap: wrap; /* Allow items to wrap to the next line if needed */
            align-items: flex-start; /* Align items at the start of the flex container */
            background-color: #F5DEB3; /* Change the background color to light tan*/
        }

        .column { /*default to stacking columns horizontally side-by-side*/
            flex: 1; /* Each column takes up equal space */
            max-width: 50%; /* Limit each column to half the width of the page */
            margin-right: 0; /* Remove the margin between columns */
            box-sizing: border-box; /* Include padding and border in the total width */
        }
        
        @media screen and (max-width: 600px) { /*when screens are <=600 pixels, stack columns vertically*/
            .column {
                width: 100%;
                float: none;
            }
        }
        
        .greeting {
            font-size: 1.2em; /* Increase the font size of the greeting */
            margin-bottom: 20px; /* Add margin to the bottom of the greeting */
        }      

        p {
            text-align: justify; /* Justify the text */
        }

        h2 {
            clear: right; /* Clear the float to ensure headings start below the image */
        }

        /* Style to remove bullets from all lists */
        ul {
            list-style-type: none;
            padding: 0; /* Remove default padding */
        }

        .columnWrapper {
            display: flex;
            flex-direction: column; /* Stack image and text vertically */
            align-items: center; /* Center items horizontally */
            text-align: center; /* Center text */
        }

        .imageContainer {
            margin-bottom: 20px; /* Add space between image and text */
        }

        img {
            display: block;
            width: 50%;
            height: auto;
            margin: 0 auto;
        }
        
    </style>
</head>
<body>

    <div class="column">
        <!-- Greeting at the top of the left column -->
        <p class="greeting">
            Hi there! My name is Elizabeth Thompson. I am a research and teaching assistant at Washington State University. I love teaching Calculus, and my research is in Topological Data Analysis and Machine Learning!
            <!-- Link to the PDF file -->
            <a href="Graduate CV.pdf" target="_blank">Here is a link to my CV.</a>
        </p>

        <!-- Education -->
        <h2>Education</h2>
        <ul>
            <li>PhD in Mathematics, Washington State University (WSU), 2021-Present</li>
            <li>BS in Mathematics & Secondary Education, Linfield University, 2017-2021</li>
        </ul>

        <!-- Presentations -->
        <h2>Presentations</h2>
        <ul>
            <li><a href="https://www.westerncriminology.org/documents/conference_proceedings/WSC_2024_Conference_Program.pdf"> Western Society of Criminology 50th Annual Conference</a> "Perceptions of Force Severity, Operationalizing Reasonableness With a Repeated Measures Design" by David Makin, Elizabeth Thompson, Mary McMillin, Bala Krishnamoorthy, and Dale Willits, Long Beach CA, February 9, 2024</li>
            <li><a href="https://www.youtube.com/watch?v=iu7vNOukW00&t=3550s">What Can Donuts Tell Us About Body Cam Data?</a> 2nd Place Three Minute Thesis Finals, WSU, March 29, 2023</li>        
            <li><a href="https://georgefoxuniversity.regfox.com/pacific-northwest-section-of-the-mathematical-association-of-america"> PNW MAA Meeting</a> "Motivating Student Learning with Mathematical Modeling in Calculus for Life Science" by Alexander Dimitrov and Elizabeth Thompson, George Fox University, March 18, 2023</li>
        </ul>

        <!-- Teaching -->
        <h2>Teaching</h2>
        <ul>
            <li>Math 202: Calculus for Business & Economics Lecture, Fall 2023</li>
            <li>Math 103: Algebra Methods & Functions Lecture, Summer 2023</li>
            <li>Math 108: Trigonometry Lecture, Summer 2022 & Summer 2023</li>
            <li>Math 140 Lab: Calculus for Life Scientists, Fall 2023-Spring 2023</li>
            <li>Math 171 Global Campus: Calculus I Lecture & Lab, Summer 2022</li>
        </ul>

        <!-- Papers -->
        <h2>Papers (in progress)</h2>
        <ul>
            <li> <a href="DQE_Final_Draft.pdf" target="_blank">Comparing TDA Methods of Time Series Analysis</a> 
                 <a href="R.Code Main.R" target="_blank">with associated R code</a> </li>
            <li> <a href="Perceptions_of_Force_Update_March_2024.pdf" target="_blank">Predictive Models for Perceptions of Police-Community Interactions Involving Use of Force</a> 
                 <a href="Perceptions_UofF_Code.ipynb" target="_blank"> with associated Python code</a> </li>
            <li> <a href="Mapper_Analysis.pdf" target="_blank">Factors That Impact Uses of Force in Police-Community Interactions</a> 
                 <a href="UofF Mapper.ipynb" target="_blank">with associated Python code </a> </li>
        </ul>
    </div>

    <div class="column">
        <!-- Image and Contact Information -->
        <div class="columnWrapper">
            <div class="imageContainer">
                <img src="https://github.com/ElizabethThompson98/ElizabethThompson98.github.io/blob/main/Directory_Photo.jpg?raw=true" alt=""/>
            </div>
            <div class="contactInfo">
                <p>
                    <strong>Contact Information</strong><br>
                    Email: elizabeth.thompson1@wsu.edu<br>
                    Office Address:<br>
                    Vancouver Undergraduate Building<br>
                    Room 251<br>
                    Washington State University<br>
                    Vancouver WA, 98661
                </p>
            </div>
        </div>
    </div>

</body>
</html>
