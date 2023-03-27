namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// Bulk1 最大SNP-indexルール
    /// </summary>
    internal class Bulk1MaxSnpIndexRule : IQtlSeqTargetVariantRule
    {
        private readonly MaxBulkSnpIndexThreshold _threshold;

        /// <summary>
        /// Bulk1 最大SNP-indexルール インスタンスを作成する。
        /// </summary>
        /// <param name="threshold">しきい値</param>
        public Bulk1MaxSnpIndexRule(MaxBulkSnpIndexThreshold threshold)
        {
            _threshold = threshold;
        }

        public bool Ok(SnpIndexVariant variant)
        {
            return variant.Bulk1.SnpIndex.Value <= _threshold.Value;
        }
    }
}
