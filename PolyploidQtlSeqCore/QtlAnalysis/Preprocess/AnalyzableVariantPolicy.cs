namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess
{
    /// <summary>
    /// 解析可能変異ポリシー
    /// </summary>
    internal class AnalyzableVariantPolicy
    {
        private readonly IAnalyzableVariantRule[] _rules;

        /// <summary>
        /// 解析可能変異ポリシーを作成する。
        /// </summary>
        public AnalyzableVariantPolicy()
        {
            _rules = new IAnalyzableVariantRule[]
            {
                new Parent1GtHomoRule(),
                new BiallelicVariantRule(),
                new Parent1NoGapVariantRule(),
                new Parent2NoGapVariantRule(),
                new Bulk1NoGapVariantRule(),
                new Bulk2NoGapVariantRule()
            };
        }

        /// <summary>
        /// 解析可能変異かどうかを判断する。
        /// 全てのルールに適合した場合のみ解析可能とする。
        /// </summary>
        /// <param name="variant">変異</param>
        /// <returns>解析可能変異ならtrue</returns>
        public bool ComplyWithAll(VcfVariant variant)
        {
            return _rules.All(x => x.Ok(variant));
        }
    }
}
