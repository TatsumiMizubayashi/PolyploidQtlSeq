using PolyploidQtlSeqCore.QtlAnalysis.IO;

namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// SnpEffアノテーションコレクション
    /// </summary>
    internal class SnpEffAnnotations
    {
        private const string DELIMITER = " | ";

        private readonly SnpEffAnnotation[] _annotations;

        /// <summary>
        /// 空のSnpEffアノテーションコレクションを作成する。
        /// </summary>
        public SnpEffAnnotations()
        {
            _annotations = Array.Empty<SnpEffAnnotation>();
            TopImpact = Impact.None;
            IsEmpty = true;
        }


        /// <summary>
        /// SnpEffアノテーションコレクションを作成する。
        /// </summary>
        /// <param name="annotations">アノテーション</param>
        public SnpEffAnnotations(SnpEffAnnotation[] annotations)
        {
            _annotations = annotations;

            TopImpact = _annotations.Length == 0
                ? Impact.None
                : _annotations.Select(x => x.Impact).Distinct().OrderBy(x => x).First();
            IsEmpty = !_annotations.Any();
        }

        /// <summary>
        /// 影響力が一番強いImpactを取得する。
        /// </summary>
        public Impact TopImpact { get; }

        /// <summary>
        /// アノテーション情報が存在するかどうかを取得する。
        /// </summary>
        public bool IsEmpty { get; }

        /// <summary>
        /// 表示アノテーションのみのSnpEffアノテーションコレクションに変換する。
        /// </summary>
        /// <param name="displayImpacts">表示するImpact</param>
        /// <returns>表示SnpEffアノテーションコレクション</returns>
        public SnpEffAnnotations ToDisplaySnpEffAnnotations(DisplayAnnotationImpacts displayImpacts)
        {
            var displayAnns = _annotations.Where(x => x.Display(displayImpacts)).ToArray();

            return displayAnns.Length == 0
                ? new SnpEffAnnotations()
                : new SnpEffAnnotations(displayAnns);
        }

        /// <summary>
        /// 出力用Impact情報に変換する。
        /// </summary>
        /// <returns>Impact情報</returns>
        public string ToImpactInfo()
        {
            var infoQuery = _annotations.Select(x => x.ToImpactInfo());

            return string.Join(DELIMITER, infoQuery);
        }

        /// <summary>
        /// 出力用Annotation情報に変換する。
        /// </summary>
        /// <returns>Annotation情報</returns>
        public string ToAnnotationInfo()
        {
            var infoQuery = _annotations.Select(x => x.ToAnnotationInfo());

            return string.Join(DELIMITER, infoQuery);
        }

        /// <summary>
        /// 出力用HGVS.C情報に変換する。
        /// </summary>
        /// <returns>HGVS.C情報</returns>
        public string ToHgvsCInfo()
        {
            var infoQuery = _annotations.Select(x => x.ToHgvsCInfo());

            return string.Join(DELIMITER, infoQuery);
        }

        /// <summary>
        /// 出力用HGVS.P情報に変換する。
        /// </summary>
        /// <returns>HGVS.P情報</returns>
        public string ToHgvsPInfo()
        {
            var infoQuery = _annotations
                .Select(x => x.ToHgvsPInfo())
                .Where(x => !string.IsNullOrEmpty(x));

            return string.Join(DELIMITER, infoQuery);
        }

    }
}
