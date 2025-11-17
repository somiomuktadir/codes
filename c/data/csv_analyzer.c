#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_RECORDS 100
#define FILENAME "student_data.csv"
#define LINE_LENGTH 256
#define TOP_BOTTOM_PERCENTAGE 5

typedef struct {
    int height;
    int weight;
    char gender;
    float CGPA;
    int age;
} StudentRecord;

StudentRecord data[MAX_RECORDS];
int recordCount = 0;
int dataLoaded = 0;

int compare_floats(const void *a, const void *b) {
    float fa = *(const float*)a;
    float fb = *(const float*)b;
    if (fa < fb) return -1;
    if (fa > fb) return 1;
    return 0;
}

int compare_records_cgpa(const void *a, const void *b) {
    StudentRecord *rec_a = (StudentRecord *)a;
    StudentRecord *rec_b = (StudentRecord *)b;
    if (rec_a->CGPA < rec_b->CGPA) return -1;
    if (rec_a->CGPA > rec_b->CGPA) return 1;
    return 0;
}

void read_csv() {
    FILE *file = fopen(FILENAME, "r");
    char line[LINE_LENGTH];
    char *token;
    int i = 0;

    if (file == NULL) {
        printf("\nERROR: File %s not found. Ensure it is in the same directory.\n", FILENAME);
        dataLoaded = 0;
        return;
    }

    if (fgets(line, LINE_LENGTH, file) == NULL) {
        fclose(file);
        printf("\nERROR: CSV file is empty.\n");
        dataLoaded = 0;
        return;
    }

    while (fgets(line, LINE_LENGTH, file) != NULL && i < MAX_RECORDS) {
        token = strtok(line, ",");
        if (token) data[i].height = atoi(token);

        token = strtok(NULL, ",");
        if (token) data[i].weight = atoi(token);

        token = strtok(NULL, ",");
        if (token) data[i].gender = token[0];

        token = strtok(NULL, ",");
        if (token) data[i].CGPA = atof(token);

        token = strtok(NULL, ",");
        if (token) data[i].age = atoi(token);

        i++;
    }

    recordCount = i;
    fclose(file);
    dataLoaded = 1;

    printf("\nSUCCESS: Read %d records from %s.\n", recordCount, FILENAME);

    if (dataLoaded) {
        printf("\n--- Loaded Data Snapshot (All %d Records) ---\n", recordCount);
        printf("+------+---------+---------+--------+------+\n");
        printf("|  ID  | HEIGHT  | WEIGHT  | GENDER | CGPA |\n");
        printf("+------+---------+---------+--------+------+\n");

        for (int j = 0; j < recordCount; j++) {
            printf("| %4d | %6d cm | %6d kg | %5c  | %4.2f |\n",
                   j + 1, data[j].height, data[j].weight, data[j].gender, data[j].CGPA);
        }

        printf("+------+---------+---------+--------+------+\n");
    }
}

void analyze_bmi() {
    int i;
    int underweight = 0;
    int normal = 0;
    int overweight = 0;
    int obese = 0;
    float bmi;
    float height_m;

    printf("\n--- 01. Body Mass Index (BMI) Report: Wellness Check ---\n");
    printf("Category | Count | Percentage\n");
    printf("---------|-------|------------\n");

    for (i = 0; i < recordCount; i++) {
        height_m = (float)data[i].height / 100.0;
        bmi = (float)data[i].weight / (height_m * height_m);

        if (bmi < 18.5) underweight++;
        else if (bmi >= 18.5 && bmi < 24.9) normal++;
        else if (bmi >= 25.0 && bmi < 29.9) overweight++;
        else obese++;
    }

    printf("Underweight| %5d | %9.2f%%\n", underweight, (float)underweight * 100.0 / recordCount);
    printf("Normal     | %5d | %9.2f%%\n", normal, (float)normal * 100.0 / recordCount);
    printf("Overweight | %5d | %9.2f%%\n", overweight, (float)overweight * 100.0 / recordCount);
    printf("Obese      | %5d | %9.2f%%\n", obese, (float)obese * 100.0 / recordCount);
    printf("---------------------------\n");
}

void analyze_quartiles() {
    float cgpa_array[MAX_RECORDS];
    int i;
    float sum = 0.0;
    float mean;
    float q1, median, q3;

    for (i = 0; i < recordCount; i++) {
        cgpa_array[i] = data[i].CGPA;
        sum += data[i].CGPA;
    }

    qsort(cgpa_array, recordCount, sizeof(float), compare_floats);

    mean = sum / recordCount;

    if (recordCount % 2 == 0) {
        median = (cgpa_array[recordCount / 2 - 1] + cgpa_array[recordCount / 2]) / 2.0;
    } else {
        median = cgpa_array[recordCount / 2];
    }

    q1 = cgpa_array[(int)(0.25 * recordCount)];
    q3 = cgpa_array[(int)(0.75 * recordCount)];

    printf("\n--- 02. Academic Performance Quartiles: Statistical Grade Snapshot ---\n");
    printf("Total Records: %d\n", recordCount);
    printf("Mean CGPA (Average): %.3f\n", mean);
    printf("Median CGPA (Q2 - 50th Pctl): %.3f\n", median);
    printf("First Quartile (Q1 - 25th Pctl): %.3f\n", q1);
    printf("Third Quartile (Q3 - 75th Pctl): %.3f\n", q3);
    printf("Interquartile Range (IQR): %.3f\n", q3 - q1);
    printf("--------------------------------------\n");
}

void analyze_correlation() {
    int i;
    float sum_h = 0.0, sum_w = 0.0;
    float mean_h, mean_w;
    float num = 0.0, den_h = 0.0, den_w = 0.0;
    float corr;

    for (i = 0; i < recordCount; i++) {
        sum_h += data[i].height;
        sum_w += data[i].weight;
    }

    mean_h = sum_h / recordCount;
    mean_w = sum_w / recordCount;

    for (i = 0; i < recordCount; i++) {
        float dev_h = data[i].height - mean_h;
        float dev_w = data[i].weight - mean_w;
        num += dev_h * dev_w;
        den_h += dev_h * dev_h;
        den_w += dev_w * dev_w;
    }

    if (den_h * den_w == 0) {
        corr = 0.0;
    } else {
        corr = num / sqrt(den_h * den_w);
    }

    printf("\n--- 03. Height vs. Weight Correlation (Pearson R): Linear Relationship Check ---\n");
    printf("Correlation Coefficient (R): %.4f\n", corr);
    printf("\nInterpretation:\n");
    if (corr > 0.7) printf("Strong positive linear relationship: Taller students tend to be heavier.\n");
    else if (corr > 0.3) printf("Moderate positive linear relationship.\n");
    else if (corr > -0.3) printf("Weak or no linear relationship.\n");
    else if (corr > -0.7) printf("Moderate negative linear relationship.\n");
    else printf("Strong negative linear relationship.\n");
    printf("--------------------------------------------------------------------------------\n");
}

void analyze_gender() {
    int i;
    int male_count = 0;
    int female_count = 0;
    float male_cgpa_sum = 0.0;
    float female_cgpa_sum = 0.0;

    for (i = 0; i < recordCount; i++) {
        if (data[i].gender == 'M' || data[i].gender == 'm') {
            male_count++;
            male_cgpa_sum += data[i].CGPA;
        } else if (data[i].gender == 'F' || data[i].gender == 'f') {
            female_count++;
            female_cgpa_sum += data[i].CGPA;
        }
    }

    printf("\n--- 04. Gender & Academic Distribution: Comparative Study ---\n");
    printf("Total Records: %d\n", recordCount);
    printf("\nMale Students:\n");
    printf("  Count: %d (%.2f%%)\n", male_count, (float)male_count * 100.0 / recordCount);
    if (male_count > 0) printf("  Average CGPA: %.3f\n", male_cgpa_sum / male_count);
    else printf("  Average CGPA: N/A\n");

    printf("\nFemale Students:\n");
    printf("  Count: %d (%.2f%%)\n", female_count, (float)female_count * 100.0 / recordCount);
    if (female_count > 0) printf("  Average CGPA: %.3f\n", female_cgpa_sum / female_count);
    else printf("  Average CGPA: N/A\n");
    printf("--------------------------------------\n");
}

void analyze_age_segmentation() {
    int i;
    int group1_count = 0, group2_count = 0, group3_count = 0;
    float group1_cgpa_sum = 0.0, group2_cgpa_sum = 0.0, group3_cgpa_sum = 0.0;

    printf("\n--- 05. Age Group CGPA Segmentation: Performance by Cohort ---\n");
    printf("Groups are: 18-20 (Young), 21-23 (Middle), 24+ (Senior)\n");

    for (i = 0; i < recordCount; i++) {
        if (data[i].age >= 18 && data[i].age <= 20) {
            group1_count++;
            group1_cgpa_sum += data[i].CGPA;
        } else if (data[i].age >= 21 && data[i].age <= 23) {
            group2_count++;
            group2_cgpa_sum += data[i].CGPA;
        } else {
            group3_count++;
            group3_cgpa_sum += data[i].CGPA;
        }
    }

    printf("\nGroup 1 (18-20 years):\n");
    printf("  Count: %d\n", group1_count);
    if (group1_count > 0) printf("  Average CGPA: %.3f\n", group1_cgpa_sum / group1_count);

    printf("\nGroup 2 (21-23 years):\n");
    printf("  Count: %d\n", group2_count);
    if (group2_count > 0) printf("  Average CGPA: %.3f\n", group2_cgpa_sum / group2_count);

    printf("\nGroup 3 (24+ years):\n");
    printf("  Count: %d\n", group3_count);
    if (group3_count > 0) printf("  Average CGPA: %.3f\n", group3_cgpa_sum / group3_count);

    printf("-------------------------------------\n");
}

// NEW ANALYSIS 6: CGPA vs. Height Correlation
void analyze_cgpa_vs_height() {
    int i;
    float sum_h = 0.0, sum_g = 0.0;
    float mean_h, mean_g;
    float num = 0.0, den_h = 0.0, den_g = 0.0;
    float corr;

    for (i = 0; i < recordCount; i++) {
        sum_h += data[i].height;
        sum_g += data[i].CGPA;
    }

    mean_h = sum_h / recordCount;
    mean_g = sum_g / recordCount;

    for (i = 0; i < recordCount; i++) {
        float dev_h = data[i].height - mean_h;
        float dev_g = data[i].CGPA - mean_g;
        num += dev_h * dev_g;
        den_h += dev_h * dev_h;
        den_g += dev_g * dev_g;
    }

    if (den_h * den_g == 0) {
        corr = 0.0;
    } else {
        corr = num / sqrt(den_h * den_g);
    }

    printf("\n--- 06. CGPA vs. Height Correlation: Is there a physical predictor for grades? ---\n");
    printf("Correlation Coefficient (R): %.4f\n", corr);
    printf("\nInterpretation:\n");
    if (fabs(corr) < 0.1) printf("Extremely weak linear relationship. Height is not a predictor for CGPA.\n");
    else if (corr > 0) printf("Weak positive relationship. Taller students have a slightly higher CGPA.\n");
    else printf("Weak negative relationship. Taller students have a slightly lower CGPA.\n");
    printf("----------------------------------------------------------------------------------\n");
}

// NEW ANALYSIS 7: Outlier Detection
void analyze_outliers() {
    int i;
    int num_outliers = (recordCount * TOP_BOTTOM_PERCENTAGE) / 100;
    StudentRecord sorted_data[MAX_RECORDS];

    for (i = 0; i < recordCount; i++) {
        sorted_data[i] = data[i];
    }
    qsort(sorted_data, recordCount, sizeof(StudentRecord), compare_records_cgpa);

    printf("\n--- 07. CGPA Outlier Report: Top & Bottom %d%% Performers ---\n", TOP_BOTTOM_PERCENTAGE);
    printf("Total records: %d. Showing Top and Bottom %d records.\n", recordCount, num_outliers);
    printf("+-------+-------+------+\n");
    printf("|  Type | CGPA  | Age  |\n");
    printf("+-------+-------+------+\n");

    printf("Top Performers:\n");
    for (i = recordCount - 1; i >= recordCount - num_outliers; i--) {
        printf("| TOP   | %5.2f | %4d |\n", sorted_data[i].CGPA, sorted_data[i].age);
    }
    printf("+-------+-------+------+\n");

    printf("Bottom Performers:\n");
    for (i = 0; i < num_outliers; i++) {
        printf("| BOTTOM| %5.2f | %4d |\n", sorted_data[i].CGPA, sorted_data[i].age);
    }
    printf("+-------+-------+------+\n");
}

// NEW ANALYSIS 8: Weight Bands
void analyze_weight_bands() {
    int i;
    int light_count = 0, medium_count = 0, heavy_count = 0;
    float light_cgpa_sum = 0.0, medium_cgpa_sum = 0.0, heavy_cgpa_sum = 0.0;
    
    // Determine bands based on overall weight distribution
    float max_w = 0, min_w = 1000;
    for (i = 0; i < recordCount; i++) {
        if (data[i].weight > max_w) max_w = data[i].weight;
        if (data[i].weight < min_w) min_w = data[i].weight;
    }
    float range = max_w - min_w;
    float band1_max = min_w + range / 3.0;
    float band2_max = min_w + 2 * range / 3.0;

    for (i = 0; i < recordCount; i++) {
        if (data[i].weight < band1_max) {
            light_count++;
            light_cgpa_sum += data[i].CGPA;
        } else if (data[i].weight < band2_max) {
            medium_count++;
            medium_cgpa_sum += data[i].CGPA;
        } else {
            heavy_count++;
            heavy_cgpa_sum += data[i].CGPA;
        }
    }

    printf("\n--- 08. Weight-Group CGPA Averages ---\n");
    printf("Weight Range | Count | Average CGPA\n");
    printf("-------------|-------|--------------\n");

    printf("Light (<%.1f kg)| %5d | %12.3f\n", band1_max, light_count, light_count > 0 ? light_cgpa_sum / light_count : 0.0);
    printf("Medium (<%.1f kg)| %5d | %12.3f\n", band2_max, medium_count, medium_count > 0 ? medium_cgpa_sum / medium_count : 0.0);
    printf("Heavy (>=%.1f kg)| %5d | %12.3f\n", band2_max, heavy_count, heavy_count > 0 ? heavy_cgpa_sum / heavy_count : 0.0);
    printf("-----------------------------------\n");
}

// NEW ANALYSIS 9: CGPA Score Distribution
void analyze_cgpa_bands() {
    int i;
    int band1 = 0; // < 3.0 (Low)
    int band2 = 0; // 3.0 - 3.49 (Mid)
    int band3 = 0; // 3.5 - 3.79 (High)
    int band4 = 0; // 3.8 - 4.0 (Top)

    for (i = 0; i < recordCount; i++) {
        if (data[i].CGPA < 3.0) band1++;
        else if (data[i].CGPA < 3.5) band2++;
        else if (data[i].CGPA < 3.8) band3++;
        else band4++;
    }

    printf("\n--- 09. CGPA Score Distribution: Grade Frequency ---\n");
    printf("CGPA Range | Status | Count | Percentage\n");
    printf("-----------|--------|-------|------------\n");
    printf(" < 3.0     | Low    | %5d | %9.2f%%\n", band1, (float)band1 * 100.0 / recordCount);
    printf(" 3.0 - 3.49| Mid    | %5d | %9.2f%%\n", band2, (float)band2 * 100.0 / recordCount);
    printf(" 3.5 - 3.79| High   | %5d | %9.2f%%\n", band3, (float)band3 * 100.0 / recordCount);
    printf(" 3.8 - 4.0 | Top    | %5d | %9.2f%%\n", band4, (float)band4 * 100.0 / recordCount);
    printf("---------------------------------------\n");
}

// NEW ANALYSIS 10: Summary Statistics
void analyze_summary_stats() {
    int i;
    float sum_h=0, sum_w=0, sum_a=0, sum_g=0;
    float min_h=1000, max_h=0, min_w=1000, max_w=0, min_a=1000, max_a=0, min_g=1000, max_g=0;

    for (i = 0; i < recordCount; i++) {
        sum_h += data[i].height;
        if (data[i].height < min_h) min_h = data[i].height;
        if (data[i].height > max_h) max_h = data[i].height;

        sum_w += data[i].weight;
        if (data[i].weight < min_w) min_w = data[i].weight;
        if (data[i].weight > max_w) max_w = data[i].weight;

        sum_a += data[i].age;
        if (data[i].age < min_a) min_a = data[i].age;
        if (data[i].age > max_a) max_a = data[i].age;

        sum_g += data[i].CGPA;
        if (data[i].CGPA < min_g) min_g = data[i].CGPA;
        if (data[i].CGPA > max_g) max_g = data[i].CGPA;
    }

    printf("\n--- 10. Summary Statistics: Range and Central Tendency ---\n");
    printf("Field | Minimum | Maximum | Range | Average\n");
    printf("------|---------|---------|-------|---------\n");
    printf("Height| %7.0f | %7.0f | %5.0f | %7.2f\n", min_h, max_h, max_h - min_h, sum_h / recordCount);
    printf("Weight| %7.0f | %7.0f | %5.0f | %7.2f\n", min_w, max_w, max_w - min_w, sum_w / recordCount);
    printf("Age   | %7.0f | %7.0f | %5.0f | %7.2f\n", min_a, max_a, max_a - min_a, sum_a / recordCount);
    printf("CGPA  | %7.2f | %7.2f | %5.2f | %7.3f\n", min_g, max_g, max_g - min_g, sum_g / recordCount);
    printf("----------------------------------------------\n");
}

// NEW ANALYSIS 11: Age by CGPA performance
void analyze_bivariate_avg_age() {
    int i;
    float total_cgpa_sum = 0.0;
    for (i = 0; i < recordCount; i++) {
        total_cgpa_sum += data[i].CGPA;
    }
    float overall_mean_cgpa = total_cgpa_sum / recordCount;

    int high_cgpa_count = 0;
    float high_cgpa_age_sum = 0.0;
    int low_cgpa_count = 0;
    float low_cgpa_age_sum = 0.0;

    for (i = 0; i < recordCount; i++) {
        if (data[i].CGPA >= overall_mean_cgpa) {
            high_cgpa_count++;
            high_cgpa_age_sum += data[i].age;
        } else {
            low_cgpa_count++;
            low_cgpa_age_sum += data[i].age;
        }
    }

    printf("\n--- 11. Average Age by Performance Group ---\n");
    printf("Overall Mean CGPA Threshold: %.3f\n", overall_mean_cgpa);
    printf("\nPerformance Group| Count | Average Age\n");
    printf("-----------------|-------|--------------\n");
    printf("High CGPA (>=Avg)| %5d | %11.2f years\n", high_cgpa_count, high_cgpa_count > 0 ? high_cgpa_age_sum / high_cgpa_count : 0.0);
    printf("Low CGPA (<Avg)  | %5d | %11.2f years\n", low_cgpa_count, low_cgpa_count > 0 ? low_cgpa_age_sum / low_cgpa_count : 0.0);
    printf("----------------------------------------\n");
}


void analyze_data_menu() {
    int choice;

    if (!dataLoaded) {
        printf("\nERROR: Data is not loaded. Please select option 1 first.\n");
        return;
    }

    do {
        printf("\n--- Data Analysis Sub-Menu (11 Tools) ---\n");
        printf("1. Body Mass Index (BMI) Report: Wellness Check\n");
        printf("2. Academic Performance Quartiles: Statistical Grade Snapshot\n");
        printf("3. Height/Weight Correlation: Linear Relationship Check\n");
        printf("4. Gender & Academic Distribution: Comparative Study\n");
        printf("5. Age Group CGPA Segmentation: Performance by Cohort\n");
        printf("-----------------------------------------------------------\n");
        printf("6. CGPA vs. Height Correlation: Is there a physical predictor?\n");
        printf("7. CGPA Outlier Report: Top & Bottom Performers\n");
        printf("8. Weight-Group CGPA Averages: Weight/Grade Comparison\n");
        printf("9. CGPA Score Distribution: Grade Frequency\n");
        printf("10. Summary Statistics: Range and Central Tendency\n");
        printf("11. Average Age by Performance Group\n");
        printf("-----------------------------------------------------------\n");
        printf("12. Back to Main Menu\n");
        printf("Enter your choice: ");
        if (scanf("%d", &choice) != 1) {
            choice = 0;
            while (getchar() != '\n');
        }

        switch (choice) {
            case 1: analyze_bmi(); break;
            case 2: analyze_quartiles(); break;
            case 3: analyze_correlation(); break;
            case 4: analyze_gender(); break;
            case 5: analyze_age_segmentation(); break;
            case 6: analyze_cgpa_vs_height(); break;
            case 7: analyze_outliers(); break;
            case 8: analyze_weight_bands(); break;
            case 9: analyze_cgpa_bands(); break;
            case 10: analyze_summary_stats(); break;
            case 11: analyze_bivariate_avg_age(); break;
            case 12: printf("\nReturning to main menu.\n"); break;
            default: printf("\nInvalid option. Try again.\n");
        }
    } while (choice != 12);
}

int main() {
    int choice;

    printf("CSV Data Analyzer Initializing.\n");

    do {
        printf("\n======= Main Menu =======\n");
        printf("1. Load & Display\n");
        printf("2. Analyze CSV Data\n");
        printf("3. Exit Program\n");
        printf("Enter your choice: ");

        if (scanf("%d", &choice) != 1) {
            choice = 0;
            while (getchar() != '\n');
        }

        switch (choice) {
            case 1:
                read_csv();
                break;
            case 2:
                analyze_data_menu();
                break;
            case 3:
                printf("\nExiting the Data Analyzer. Goodbye.\n");
                break;
            default:
                printf("\nInvalid option. Try again.\n");
        }
    } while (choice != 3);

    return 0;
}