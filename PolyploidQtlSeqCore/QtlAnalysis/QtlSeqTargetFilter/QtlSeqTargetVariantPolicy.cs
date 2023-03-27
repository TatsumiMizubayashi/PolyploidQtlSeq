namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// QTL-seq解析対象変異ポリシー
    /// </summary>
    internal class QtlSeqTargetVariantPolicy
    {
        private readonly IQtlSeqTargetVariantRule[] _rules;

        /// <summary>
        /// QTL-seq解析対象変異ポリシーを作成する。
        /// </summary>
        /// <param name="rules">ルール</param>
        public QtlSeqTargetVariantPolicy(IQtlSeqTargetVariantRule[] rules)
        {
            _rules = rules;
        }

        /// <summary>
        /// 解析対象変異かどうかを判断する。
        /// ルールを全て満たす変異のみが解析対象となる。
        /// </summary>
        /// <param name="variant">変異</param>
        /// <returns>解析対象変異ならtrue</returns>
        public bool ComplyWithAll(SnpIndexVariant variant)
        {
            return _rules.All(x => x.Ok(variant));
        }
    }
}
