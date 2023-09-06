namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// QTL-seq解析対象変異ポリシー設定
    /// </summary>
    internal class QtlSeqTargetPolicySettings
    {
        /// <summary>
        /// QTL-seq解析対象変異ポリシー設定を作成する。
        /// </summary>
        /// <param name="settingValue">設定値</param>
        public QtlSeqTargetPolicySettings(IQtlSeqTargetPolicySettingValue settingValue)
        {
            Parent1MostAlleleRateThreshold = new Parent1MostAlleleRateThreshold(settingValue.Parent1MostAlleleRateThreshold);
            Parent2SnpIndexRange = new Parent2SnpIndexRange(settingValue.Parent2SnpIndexRange);
            MinimumDepthThreshold = new MinimumDepthThreshold(settingValue.MinimumDepthThreshold);
            MaxBulkSnpIndexThreshold = new MaxBulkSnpIndexThreshold(settingValue.MaxBulkSnpIndexThreshold);
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

    }
}
