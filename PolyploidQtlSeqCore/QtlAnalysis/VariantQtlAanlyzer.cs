using PolyploidQtlSeqCore.QtlAnalysis.Distribution;
using PolyploidQtlSeqCore.Share;

namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// 変異のQTL解析
    /// </summary>
    internal class VariantQtlAanlyzer
    {
        private const double P95_SIGNIFICANT_LEVEL = 0.95;
        private const double P99_SIGNIFICANT_LEVEL = 0.99;

        private readonly F1NoQtlDeltaSnpIndexDistributionGenerator _distributionGenerator;
        private readonly ThreadNumber _threadNumber;

        /// <summary>
        /// 変異QTLアナライザーを作成する。
        /// </summary>
        /// <param name="distributionOption">分布オプション</param>
        /// <param name="threadNumber">スレッド数</param>
        public VariantQtlAanlyzer(NoQtlDistributionSettings distributionOption, ThreadNumber threadNumber)
        {
            _distributionGenerator = new F1NoQtlDeltaSnpIndexDistributionGenerator(distributionOption);
            _threadNumber = threadNumber;
        }

        /// <summary>
        /// 変異のQTL解析を行う。
        /// </summary>
        /// <param name="variants">SNP-index変異</param>
        /// <returns>QTL解析済み変異</returns>
        public SnpIndexVariantWithQtl[] Analyze(SnpIndexVariant[] variants)
        {
            var pOption = new ParallelOptions()
            {
                MaxDegreeOfParallelism = _threadNumber.Value
            };

            var qtlVariants = new SnpIndexVariantWithQtl[variants.Length];
            Parallel.For(0, variants.Length, pOption, i =>
            {
                var variant = variants[i];
                var noQtlDistribution = _distributionGenerator.Generate(variant.Bulk1, variant.Bulk2);
                var p95Threshold = noQtlDistribution.GetThreshold(P95_SIGNIFICANT_LEVEL);
                var p99Threshold = noQtlDistribution.GetThreshold(P99_SIGNIFICANT_LEVEL);
                var pValue = noQtlDistribution.GetPValue(variant.DeltaSnpIndex);

                qtlVariants[i] = new SnpIndexVariantWithQtl(variant, p95Threshold, p99Threshold, pValue);
            });

            return qtlVariants;
        }

    }
}
