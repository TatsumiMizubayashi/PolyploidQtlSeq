namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// QTL-seq解析対象変異ポリシー設定
    /// </summary>
    internal class QtlSeqTargetPolicySettings
    {
        private static readonly IReadOnlyDictionary<string, string> _toLongNameDictionary = new Dictionary<string, string>()
        {
            [Parent1MostAlleleRateThreshold.SHORT_NAME] = Parent1MostAlleleRateThreshold.LONG_NAME,
            [Parent1MostAlleleRateThreshold.LONG_NAME] = Parent1MostAlleleRateThreshold.LONG_NAME,

            [Parent2SnpIndexRange.SHORT_NAME] = Parent2SnpIndexRange.LONG_NAME,
            [Parent2SnpIndexRange.LONG_NAME] = Parent2SnpIndexRange.LONG_NAME,

            [MinimumDepthThreshold.SHORT_NAME] = MinimumDepthThreshold.LONG_NAME,
            [MinimumDepthThreshold.LONG_NAME] = MinimumDepthThreshold.LONG_NAME,

            [MaxBulkSnpIndexThreshold.SHORT_NAME] = MaxBulkSnpIndexThreshold.LONG_NAME,
            [MaxBulkSnpIndexThreshold.LONG_NAME] = MaxBulkSnpIndexThreshold.LONG_NAME
        };

        /// <summary>
        /// LongName変換項目を辞書に追加する。
        /// </summary>
        /// <param name="dictionary">LongName変換辞書</param>
        public static void AddLongNameKeyValuePair(Dictionary<string, string> dictionary)
        {
            foreach (var keyValuePair in _toLongNameDictionary)
            {
                dictionary.Add(keyValuePair.Key, keyValuePair.Value);
            }
        }


        /// <summary>
        /// QTL-seq解析対象変異ポリシー設定を作成する。
        /// </summary>
        /// <param name="settingValue">オプションの値</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public QtlSeqTargetPolicySettings(IQtlSeqTargetPolicySettingValue settingValue, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Parent1MostAlleleRateThreshold = new Parent1MostAlleleRateThreshold(settingValue.Parent1MostAlleleRateThreshold, parameterDictionary, userOptionDictionary);
            Parent2SnpIndexRange = new Parent2SnpIndexRange(settingValue.Parent2SnpIndexRange, parameterDictionary, userOptionDictionary);
            MinimumDepthThreshold = new MinimumDepthThreshold(settingValue.MinimumDepthThreshold, parameterDictionary, userOptionDictionary);
            MaxBulkSnpIndexThreshold = new MaxBulkSnpIndexThreshold(settingValue.MaxBulkSnpIndexThreshold, parameterDictionary, userOptionDictionary);
        }

        /// <summary>
        /// 親1 最多アレル割合のしきい値を取得する。
        /// </summary>
        public Parent1MostAlleleRateThreshold Parent1MostAlleleRateThreshold { get; }

        /// <summary>
        /// 親2 SNP-index範囲を取得する。
        /// </summary>
        public Parent2SnpIndexRange Parent2SnpIndexRange { get; }

        /// <summary>
        /// 最低Depthしきい値を取得する。
        /// </summary>
        public MinimumDepthThreshold MinimumDepthThreshold { get; }

        /// <summary>
        /// 最大Bulk SNP-indexしきい値を取得する。
        /// </summary>
        public MaxBulkSnpIndexThreshold MaxBulkSnpIndexThreshold { get; }

        /// <summary>
        /// QTLseq解析対象変異ポリシーを作成する。
        /// </summary>
        /// <returns>QTLseq解析対象変異ポリシー</returns>
        public QtlSeqTargetVariantPolicy CreatePolicy()
        {
            var rules = new IQtlSeqTargetVariantRule[]
            {
                new Parent1MostAlleleRateRule(Parent1MostAlleleRateThreshold),
                new Parent2SnpIndexRule(Parent2SnpIndexRange),
                new Parent1DepthRule(MinimumDepthThreshold),
                new Parent2DepthRule(MinimumDepthThreshold),
                new Bulk1DepthRule(MinimumDepthThreshold),
                new Bulk2DepthRule(MinimumDepthThreshold),
                new Bulk1MaxSnpIndexRule(MaxBulkSnpIndexThreshold),
                new Bulk2MaxSnpIndexRule(MaxBulkSnpIndexThreshold)
            };

            return new QtlSeqTargetVariantPolicy(rules);
        }

        /// <summary>
        /// パラメータファイルに記載する行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        public string[] ToParameterFileLines()
        {
            return new[]
            {
                Parent1MostAlleleRateThreshold.ToParameterFileLine(),
                Parent2SnpIndexRange.ToParameterFileLine(),
                MinimumDepthThreshold.ToParameterFileLine(),
                MaxBulkSnpIndexThreshold.ToParameterFileLine()
            };
        }
    }
}
