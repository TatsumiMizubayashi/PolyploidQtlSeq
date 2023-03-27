namespace PolyploidQtlSeqCore.QtlAnalysis.IO
{
    /// <summary>
    /// SNP-indexファイルクリエーター
    /// </summary>
    internal static class SnpIndexFileCreator
    {
        private const string FILE_NAME = "SNP-Index.txt";
        internal const string QTL = "QTL";

        private static readonly string[] _fieldNames = new[]
        {
            "Chr",
            "Position",
            "Ref Allele",
            "P1 Allele",
            "Bulk1 Allele",
            "Bulk2 Allele",
            "Window P99",
            "Window P95",
            "P99",
            "P95",
            "Bulk1 Depth",
            "Bulk2 Depth",
            "P99 Threshold",
            "P95 Threshold",
            "Bulk1 SNP-index",
            "Bulk2 SNP-index",
            "Delta SNP-index",
            "P-Value",
            " -log10(P)",   // Excelで開いた際に数式と認識されないようにするためスペース入れてる
            "Window P-Value",
            "Window -log10(P)",
        };

        private static readonly string[] _snpEffFieldNames = new[]
        {
            "Impact",
            "Annotation",
            "HGVS.c",
            "HGVS.p"
        };

        /// <summary>
        /// SNP-indexファイルを作成する。
        /// </summary>
        /// <param name="outDir">出力ディレクトリ</param>
        /// <param name="variants">変異情報</param>
        /// <param name="displayImpacts">表示Impact</param>
        public static void Create(OutputDirectory outDir, SnpIndexVariantWithSlidingWindowQtl[] variants,
            DisplayAnnotationImpacts displayImpacts)
        {
            var updateDisplayImpacts = UpdateDisplayAnnotationImpacts(variants, displayImpacts);

            var filePath = outDir.CreateFilePath(FILE_NAME);
            using var writer = new StreamWriter(filePath);
            var fieldLine = string.Join("\t", GetFieldNames(updateDisplayImpacts));
            writer.WriteLine(fieldLine);

            foreach (var variant in variants)
            {
                var values = GetVariantValues(variant, updateDisplayImpacts);
                var line = string.Join("\t", values);
                writer.WriteLine(line);
            }
        }


        /*
         * QTL解析のみの場合はSnpEffオプションがないためオプションの値から推定できない。
         * そのため変異情報から推定する。
         */

        /// <summary>
        /// DisplayAnntatonImapctを更新する。
        /// 変異情報にアノテーション情報が無い場合はNoneに更新する。
        /// </summary>
        /// <param name="variants">変異情報</param>
        /// <param name="displayImpacts">ユーザー設定DisplayImpacts</param>
        /// <returns>更新されたDisplayImpacts</returns>
        private static DisplayAnnotationImpacts UpdateDisplayAnnotationImpacts(SnpIndexVariantWithSlidingWindowQtl[] variants,
            DisplayAnnotationImpacts displayImpacts)
        {
            // SnpEffアノテーションは情報付加した場合は必ず何らかの情報がつけられるのでEmptyにならない
            return variants[0].Annotations.IsEmpty
                ? new DisplayAnnotationImpacts("None", new Dictionary<string, string>(), new Dictionary<string, bool>())
                : displayImpacts;
        }


        /// <summary>
        /// フィールド名を取得する。
        /// </summary>
        /// <param name="displayImpacts">表示Impact</param>
        /// <returns>フィールド名</returns>
        private static string[] GetFieldNames(DisplayAnnotationImpacts displayImpacts)
        {
            return displayImpacts.DisplayFlag == Impact.None
                ? _fieldNames
                : _fieldNames.Concat(_snpEffFieldNames).ToArray();
        }

        /// <summary>
        /// 変異の値を取得する。
        /// </summary>
        /// <param name="variant">変異</param>
        /// <param name="displayImpacts">表示Impact</param>
        /// <returns>変異の値</returns>
        private static string[] GetVariantValues(SnpIndexVariantWithSlidingWindowQtl variant, DisplayAnnotationImpacts displayImpacts)
        {
            var genomePos = variant.GenomePosition;
            var bulk1 = variant.Bulk1;
            var bulk2 = variant.Bulk2;
            var maxScoreWindowQtl = variant.MaxScoreSlidingWindowQtl;
            var annValues = GetAnnotationValues(variant, displayImpacts);

            return new[]
            {
                genomePos.ChrName,
                genomePos.Start.ToString(),
                variant.RefAllele,
                variant.Parent1.Allele,
                bulk1.Allele,
                bulk2.Allele,
                maxScoreWindowQtl.P99QtlMark(),
                maxScoreWindowQtl.P95QtlMark(),
                variant.P99Qtl ? QTL : "",
                variant.P95Qtl ? QTL : "",
                bulk1.Depth.ToString(),
                bulk2.Depth.ToString(),
                variant.P99ThresholdDeltaSnpIndex.Value.ToString(),
                variant.P95ThresholdDeltaSnpIndex.Value.ToString(),
                bulk1.SnpIndex.Value.ToString(),
                bulk2.SnpIndex.Value.ToString(),
                variant.DeltaSnpIndex.Value.ToString(),
                variant.PValue.Value.ToString(),
                variant.Score.Value.ToString(),
                maxScoreWindowQtl.PValue.ToString(),
                maxScoreWindowQtl.Score.ToString()
            }
            .Concat(annValues)
            .ToArray();
        }

        /// <summary>
        /// SnpEffアノテーション値を取得する。
        /// </summary>
        /// <param name="variant">変異</param>
        /// <param name="displayImpacts">表示Impact</param>
        /// <returns>SnpEffアノテーション値</returns>
        private static string[] GetAnnotationValues(SnpIndexVariantWithSlidingWindowQtl variant, DisplayAnnotationImpacts displayImpacts)
        {
            if (displayImpacts.DisplayFlag == Impact.None) return Array.Empty<string>();

            var displayAnnotations = variant.Annotations.ToDisplaySnpEffAnnotations(displayImpacts);
            if (displayAnnotations.IsEmpty) return Array.Empty<string>();

            return new[]
            {
                displayAnnotations.ToImpactInfo(),
                displayAnnotations.ToAnnotationInfo(),
                displayAnnotations.ToHgvsCInfo(),
                displayAnnotations.ToHgvsPInfo()
            };
        }
    }
}
